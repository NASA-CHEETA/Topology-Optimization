import numpy as np    


def get_passive_element_list(file, output_file, skiplines, writeFlag):
    '''
    Reads the list of elements  (boundary, surface, etc),
    parses them to remove commas and spaces, and writes to a file if required.
    Mainly designed to read elements to assign passive elements, but can be used for 
    other purposes as well
    '''

    print('File loading ...')
    with open(file,"r") as f:
        data = f.readlines()


    c=[]
    count = 0

    print('parsing data ...')
    for linenum, line in enumerate(data[skiplines:]):
        b = line.split()
        for each in b:
            c.append(int(each.replace(',','')))
            #print(c[count])
            count+=1
    
    print('Writing file ...')
    if writeFlag:
        np.savetxt(output_file, c, fmt='%i')
        print('Complete!')

    return 0

if __name__ == "__main__":

    file = 'sparelements.nam'    
    output_file = 'sparElementList.nam'
    skiplines= 1    # Line starts at 1
    writeFlag = True

    c = get_passive_element_list(file, output_file, skiplines, writeFlag)
    print()
    assert (c == 0)