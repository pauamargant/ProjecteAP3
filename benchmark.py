#####################
#   Pau Amargant    #
#####################

#
# Program intended to test the AP3 Project
#
# USAGE: python3 bench.py executable_file
# ex: python3 bench.py exh.exe



import subprocess
import sys



#IMPORTANT VARIABLES
public_benchs: str = "public_benchs/" #Public benchmarks directory
temp_folder: str = "temp" #temporary folder
  
max_time = 60 #max time for each test


difs=["easy","med","hard"] #names of the difficulties
n_dif = {"easy":10,"med":10, "hard":20} #number of test of each difficulty


prev_res = { 'easy': [14, 0, 15, 0, 0, 5, 0, 3, 10, 1], 
'med': [10, 15, 39, 27, 17, 28, 14, 12, 12, 10], 
'hard': [13, 0, 0, 8, 66, 0, 0, 0, 140, 0, 7, 4, 34, 15, 1, 109, 11, 1, 0, 29]}

if __name__ == "__main__":
    program: str = "./"+sys.argv[1]
    output: str = "temp_sol.txt"

    resultats = {"easy": [], "med":[],"hard":[]} 
    pens = {"easy": [], "med":[],"hard":[]} 
    print("STARTING TEST for "+program+" with timeout at "+str(max_time))
    #for each difficulty we test all tests
    for dif in difs:
        for i in range(n_dif[dif]):
            print(f"#\n# {dif} {i+1}\n#\n")
            index = i+1
            #we execute the program for max_time seconds
            try:
                subprocess.run([program,public_benchs+dif+"-"+str(index)+".txt",output],timeout=max_time)
            except subprocess.TimeoutExpired:
                print(f'Timeout ({max_time}s) expired')
            subprocess.run(["./check.exe",public_benchs+dif+"-"+str(index)+".txt",output])
            current_pen = "-1"
            #we check the obtained result and store it
            with open(output) as f:
                lines = f.readlines()
                current_pen = str(lines[0].split()[0])
                current_time = str(lines[0].split()[1])
            print("Penalitzacio obtinguda es "+ current_pen)
            resultats[dif].append((f" t={current_time} {current_pen}",prev_res[dif][i]))
            pens[dif].append(int(current_pen))
            print("\n")
    print("\n\n IMPRIMIM RESULTAT:")
    print("|index|tested difficulty|stored result|")
    for dif in difs:
        print(f"|   |{dif}|   |")
        for i in range(len(resultats[dif])):
            print(f"| {i+1} | {resultats[dif][i][0]} | {resultats[dif][i][1]} |")
    print(pens)



    

