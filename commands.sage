load("orientations.sage")

def extract_theta_chars_from_list_with_counting(l):
    for i in range(len(l)):
        if i % 50 == 0:
            print(i)
        theta_char = l[i].get_theta_char_orientation()
        if not l[i].is_theta_char(theta_char):
            print("error at " + str(i))
            print("g6 string ") + l[i].graph6_string()

def graph_list(n):
    return [CycleCocycleSystem(G) for G in graphs(n) if G.is_biconnected()]
