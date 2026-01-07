"""
Questions Based on Algorithms by Dasgupta, Papadimitriou, Vazirani. McGraw-Hill. 2006
"""

# FIBONACCI NUMBERS
"""
Calculating the Fibonacci Number
:param fname: value, n, that's less than or equal to 25
:results: fibonacci number (F_n)  
"""
def fib_num(n): 
    if n == 0:
        return 0
    elif n == 1: 
        return 1
    else: 
        return fib_num(n-1) + fib_num(n-2)
 
    
# BINARY SEARCH 
"""
Binary Search
:param fname: file containing all information needed for the initial call of the aux function
    :n: positive integer <= 10^5 (rep len of the sorted array)
    :m: positive integer <= 10^5 (rep len of the list)
    :sorted_arr: A[1,...,n] range(-10^5, 10^5 + 1) of integers
    :lst: a list of m integers (-10^5 <= k1, k2,..., km <=10^5)
:results: index for each ki from the list of m integers as found in the sorted array
"""
def binary_search(fname): 
    with open(fname, "r") as file: 
        lines = file.read().splitlines()
        # information not needed
        n = int(lines[0])
        # m = int(lines[1])
        sorted_arr = list(map(int, lines[2].split()))
        lst = list(map(int, lines[3].split()))
        
        indices_lst = []
        for target in lst: 
            index = binary_search_aux(sorted_arr, target, 0, n)
            indices_lst.append(index)
        return indices_lst
    
"""
binary_search_aux (helper function)
recursively search for the target value in the sorted array
:param sorted_arr: A[1,...,n] range(-10^5, 10^5 + 1) of integers
:param target: integer that is being searched for
:results: index of the target in the sorted array or -1 if its not found
"""     
def binary_search_aux(arr, target, low, high): 
    if low > high: 
        return -1
    # calculate the middle of the sorted array
    mid = (low + high) // 2
    if arr[mid] == target: 
        return mid + 1
    elif arr[mid] > target: 
        return binary_search_aux(arr, target, low, mid-1)
    else: 
        return binary_search_aux(arr, target, mid+1, high)            


# DEGREE ARRAY
"""

"""
def degree_array(fname): 
    return

if __name__ == "__main__":
    # fibonacci number 
    fib_idx = int(open("rosalind_fibo.txt", "r").read())
    print(fib_num(fib_idx))
    
    # binary search 
    indices_lst = binary_search("rosalind_bins.txt")
    # convert result list int values to string values
    # combine values with space seperator and write to file
    with open("ros_bins_result.txt", "w") as file: 
        file.write(" ".join(map(str, indices_lst)))