

def get_average_runtime(benchmark_method, iterations, dataset):
    
    temp = "somemethod"

    for _ in range(iterations):

        #TEMP IS TEMPORARY! FIX LOGIC FOR PASSING METHOD AS ARGUMENT LATER!!!!!!!!!

        interval = temp("data/benchmark_sequences/astral_10000.txt") 
        sum += interval


    return sum / iterations

def get_total_dataset_residues(dataset_sequences):

    temp_s = 0

    for seq in dataset_sequences:
        temp_s += len(seq)
    
    return temp_s
