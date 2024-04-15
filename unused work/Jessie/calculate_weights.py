def calculate_weights(class_freq, total_samples):
    # METHOD 1
    weights_1 = []
    for value in class_freq.values():
        weight = total_samples/value       
        weights_1.append(weight)
    print("Method 1:", weights_1)

    # METHOD 2
    weights_2 = []
    for value in class_freq.values():
        weight = total_samples/(value*len(class_freq))
        weights_2.append(weight)
    print("Method 2:", weights_2)

    # METHOD 3
    weights_3 = []
    for value in class_freq.values():
        weight = max(class_freq.values())/value
        weights_3.append(weight)
    print("Method 3:", weights_3)

    # METHOD 4 (final choice method)
    weights_4 = []
    for value in class_freq.values():
        weight = (1/value)
        weights_4.append(weight)
    print("Method 4:", weights_4)                        



def txt_to_data(file):
    class_freq = {}
    with open(file) as data:
        for line in data:
            k, v = line.strip().split('  ')
            class_freq[k] = int(v)
    print(class_freq)
    total_samples = sum(class_freq.values())
    calculate_weights(class_freq, total_samples)


txt_to_data('res_dict.txt')


