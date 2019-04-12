import sys
import time
from multiprocessing import Pool
import multiprocessing
import tqdm
import matplotlib.pyplot as plt

# check for correct number of arguments
if (len(sys.argv) < 2):
	print("Too few arguments")
	sys.exit()

# declare arr and rRNA_arr as global variables
global arr
global rRNA_arr

# set arr and rRNA_arr as empty arrays
rRNA_arr = []
arr = []

try:
	# open rRNA file	
	rRNA = open("EColi_rRNA.txt", "r")
	
	# write file to an array and close the file
	for line in rRNA:
		rRNA_arr.append(line)
	rRNA.close()
	
	# open file specified by command line	
	RNAseq = open(sys.argv[1], "r")
	
	# write file to an array and close the file (enumerate to allow for only parsing part of file)
	print("generating array")
	for num, line in enumerate(RNAseq, 1):
		arr.append(line)
		if(num >= (130000*4)):
			break
	RNAseq.close()
	
# print any exceptions that occur and exit	
except Exception as e:
	print(e)
	sys.exit()

# write array to new file
def write_to_file(array, file):
	with open(file, 'w+') as output:
		for i in array:
			output.write(i)
	output.close()

write_to_file(arr, 'input.fastq')

# align sequences using biopython's pairwise2 package
def align_seq(sample):
	from Bio import pairwise2
	pbar.update(multiprocessing.cpu_count()*3)
	pid = []
	for i in range(3):
		alignments = pairwise2.align.localms(rRNA_arr[i], arr[sample], 1, -3, -5, -4, one_alignment_only = True)
		id = 0
		for j in range(len(alignments[0][0])):
			if(alignments[0][0][j] == alignments[0][1][j] and alignments[0][0][j] != '-'):
				id = id + 1
		if(id >= 85):
			return(id)
		pid.append(int((id / min(len(rRNA_arr[i]), len(arr[sample]))) * 100))
	output = max(pid)
	return(output)

# create an array for storing indices of DNA from input file
input_arr = []
for i in range(1, len(arr), 4):
	input_arr.append(i)
# create a progress bar
pbar = tqdm.tqdm(total = len(input_arr) * 3)

# run samples through align_seq using pool to multithread
if __name__ == '__main__':
	with Pool() as pool:
		# get start time
		start = time.time()
		result = pool.map(align_seq, input_arr)
		count = 0
		for i in range(len(result)-1, -1, -1):
			if(result[i] >= 85):
				count = count + 1
				index = i*4
				arr.pop(index+3)
				arr.pop(index+2)
				arr.pop(index+1)
				arr.pop(index)
		print()
		print(count)
		write_to_file(arr, 'output.fastq')
		# get end time and print time taken
		end = time.time()
		print()
		print("Time taken =", int(end-start), "seconds")
		
		plt.hist(result)
		plt.xlabel("Percent Identity")
		plt.ylabel("Number of Values")
		plt.title("Max Percent Identity of SRR1020876 to E.Coli rRNA")
		plt.xlim(0, 100)
		plt.grid(True)
		plt.show()