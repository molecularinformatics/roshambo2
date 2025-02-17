#!/usr/bin/env python
#
# MIT License
# 
# Copyright (c) 2025 molecularinformatics  
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import multiprocessing
import base64
import uuid
import logging
from flask import Flask, request, jsonify
import pickle
import yaml
from roshambo2 import Roshambo2ServerMode

# Load and validate configuration
def load_and_validate_config(yaml_file):
    with open(yaml_file, 'r') as file:
        config = yaml.safe_load(file)

    # Required keys
    required_keys = ['dataset_files', 'hostname', 'port', 'api_name', 'verbosity']
    
    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        raise ValueError(f"Missing required configuration keys: {', '.join(missing_keys)}")

    return config

# Worker function to process search tasks
def worker(searcher, task_queue, task_status, shared_progress, logger):
    logger.debug("Worker running")
    while True:
        task_id, query_data, options, get_structures = task_queue.get()
        logger.debug("worker got a task")

        try:
            logger.debug(f"processing task {task_id}")
            task_status[task_id] = {'status': 'running', 'result': None, 'mols': None}
            if 'backend' not in options:
                options['backend'] = 'cuda'
            results, rdkit_binary_mols_dict = searcher.search(query_data, options, shared_progress, get_structures)

            task_status[task_id] = {'status': 'completed', 'result': results, 'mols': rdkit_binary_mols_dict}
            logger.debug(f"completed task {task_id}")
        except Exception as e:
            error_message = f"Error: {str(e)}"
            task_status[task_id] = {'status': 'failed', 'result': [error_message]}
            logger.warning(f"failed task {task_id} with error {error_message}")

# Flask app initialization and routes
def create_app(config, searcher, task_queue, task_status, shared_progress, logger):
    app = Flask(__name__)

    @app.route(config['api_name'], methods=['POST'])
    def search_dataset():
        data = request.json
        encoded_query_data = data.get('query_data')
        options = data.get('options', {})
        get_structures = data.get('get_structures', True)

        try:
            query_data = pickle.loads(base64.b64decode(encoded_query_data))
        except Exception as e:
            return jsonify({"error": f"Failed to unpickle data: {str(e)}"}), 400

        task_id = str(uuid.uuid4())
        task_status[task_id] = {'status': 'pending', 'result': None, 'mols': None}
        task_queue.put((task_id, query_data, options, get_structures))
        return jsonify({"task_id": task_id}), 202

    @app.route(config['api_name'] + '/status/<task_id>', methods=['GET'])
    def check_status(task_id):
        if task_id not in task_status:
            return jsonify({"error": "Invalid task ID"}), 404

        status_info = task_status[task_id]
        if "running" in status_info['status']:
            status_info['status'] = f"running: {int(shared_progress.value * 100)}%"
        return jsonify({"task_id": task_id, "status": status_info['status']}), 200

    @app.route(config['api_name'] + '/result/<task_id>', methods=['GET'])
    def get_result(task_id):
        if task_id not in task_status:
            return jsonify({"error": "Invalid task ID"}), 404

        status_info = task_status[task_id]
        if status_info['status'] == 'failed':
            return jsonify(status_info['result'])
        if status_info['status'] != 'completed':
            return jsonify({"error": "Task not completed yet"}), 400

        result = status_info['result']
        json_dataframes_dict = {key: df.to_json() for key, df in result.items()}
        del task_status[task_id]
        return jsonify(json_dataframes_dict), 200

    @app.route(config['api_name'] + '/results_with_structures/<task_id>', methods=['GET'])
    def get_results_and_structures(task_id):
        if task_id not in task_status:
            return jsonify({"error": "Invalid task ID"}), 404

        status_info = task_status[task_id]
        if status_info['status'] == 'failed':
            return jsonify(status_info['result'])
        if status_info['status'] != 'completed':
            return jsonify({"error": "Task not completed yet"}), 400

        result = status_info['result']
        json_dataframes_dict = {key: df.to_json() for key, df in result.items()}
        mols_dict = status_info['mols']
        base64_encoded_mols_dict = {key: [base64.b64encode(mol).decode('utf-8') for mol in mols] for key, mols in mols_dict.items()}
        data = {"molecules": base64_encoded_mols_dict, "scores": json_dataframes_dict}
        del task_status[task_id]
        return jsonify(data), 200

    return app

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python server.py <config_file.yaml>")
        sys.exit(1)

    try:
        config_file = sys.argv[1]
        config = load_and_validate_config(config_file)

        dataset_files = config['dataset_files']
        verbosity = config['verbosity']

        # Initialize searcher
        searcher = Roshambo2ServerMode(dataset_files, color=True, verbosity=verbosity, n_gpus=1)

        # Multiprocessing setup
        manager = multiprocessing.Manager()
        shared_dict = manager.dict()
        shared_queue = manager.Queue()
        shared_progress = manager.Value('f', 0.0)

        # Logger setup
        logger = logging.getLogger('roshambo2_server')
        logger.setLevel(logging.DEBUG)

        # Start worker process
        process = multiprocessing.Process(target=worker, args=(searcher, shared_queue, shared_dict, shared_progress, logger))
        process.daemon = True
        process.start()

        # Create and run the app
        app = create_app(config, searcher, shared_queue, shared_dict, shared_progress, logger)
        app.run(host=config['hostname'], port=config['port'], debug=False)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
