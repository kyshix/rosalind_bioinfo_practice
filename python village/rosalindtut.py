# TASK #2 - VARIABLES AND SOME ARITHMETIC
# lackluster method
a = 877 
b = 909
c = a**2 + b**2
print(c)

# alternatively to take in the files as god intended
"""
Square Hypotenuse of Two Values (Arithmetic Practice)

:param fname: file containing two numbers
:returns: integer corresponding to the square of the hypothenuse of a 
    right triange with lengths of the two numbers (a,b)
"""
def square_hypotenuse(fname):
    values = open(fname, "r").read().split()
    a = int(values[0])
    b = int(values[1])
    return (a**2) + (b**2)

# TASK #3 - STRINGS AND LISTS
"""
Splice Phrases by Indices
:param fname: file containing a string of length of at most 200 characters
    and 4 integers (a, b, c, d) that will be used to splice two words out
:returns: spliced phrase 
"""
def splice(fname): 
    lines =  open(fname, "r").readlines()
    text = lines[0]
    values = lines[1].split()
    # note that 1 is added to the ending indices to make it inclusive
    a = int(values[0])
    b = int(values[1]) + 1
    c = int(values[2])
    d = int(values[3]) + 1
    return text[a:b] + " " + text[c:d]

# TASK 4 - CONDITIONS AND LOOPS
"""
Sum of All Numbers From Two Values
:param fname: file containing two positive numbers (a<b<10000)
:returns: integer of the sum of values from a to b (inclusive for both)
"""
def sum_btwn(fname): 
    values = open(fname, "r").read().split()
    a = int(values[0])
    b = int(values[1])
    sum = 0
    for val in range(a,b+1): 
        if val%2 == 1: 
            sum += val
    return sum

# TASK 5 - WORKING WITH FILES
"""
Even Lines of a File
:param fname: file with at most 1000 lines (1-based number of lines)
:returns: new file with only the even lines
"""
def even_lines(fname):
    try: 
        new_fname = "tut5result.txt" 
        with open(fname, "r") as og_file: 
            with open(new_fname, "w") as result:
                count = 1
                for line in og_file:
                    if count%2 == 0:
                        result.write(line)
                    count += 1
            result.close()
        og_file.close()
                        
    except FileNotFoundError: 
        print(f"Error: the file '{fname}' was not found")
    except Exception as e:
        print(f"An error occurred: {e}")

# TASK 6 - DICTIONARIES
"""
Count of Unique Word Occurences
:param fname: file with a string of length of at most 10000 characters
:returns: number of occurences each word (case sensitive)
"""
def count_words(fname): 
    with open(fname, "r") as s:
        words = s.read().split()
    
        s_dict = {}
        for word in words: 
            if word in s_dict:
                s_dict[word] += 1
            else:
                s_dict[word] = 1
        
        tut6fname = "tut6result.txt"
        with open(tut6fname, "w") as result:
            for key, value in s_dict.items(): 
                result.write(f"{key} {value}\n")
                
if __name__ == "__main__":
    print(square_hypotenuse("rosalind_ini2.txt"))
    print(splice("rosalind_ini3.txt"))
    print(sum_btwn("rosalind_ini4.txt"))
    even_lines("rosalind_ini5.txt")
    count_words("rosalind_ini6.txt")