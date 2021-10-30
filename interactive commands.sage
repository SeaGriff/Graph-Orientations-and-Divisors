def extract_theta_chars_from_list_with_counting(l):
    for i in range(len(l)):
        div = l[i].big_theta_char_divisor()
        if i % 150 == 0:
            print(i)
        if div == False:
            print("thetaless at " + str(i))
