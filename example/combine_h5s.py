from roshambo2.prepare import combine_datasets



if __name__ == "__main__":

    combine_datasets(['dataset_split_0.h5', 'dataset_split_1.h5', 'dataset_split_2.h5'], 'combined.h5')


