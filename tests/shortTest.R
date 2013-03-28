#this is a short test
library(sciClone)

#read in vaf data
v = read.table("data/vafs.dat",header=T);
v1 = v[1:100,c(1,2,8,9,10)]
v2 = v[1:100,c(1,2,11,12,13)]
v3 = v[1:100,c(1,2,14,15,16)]

#read in regions to exclude
regions = read.table("data/exclude.chr1")

#read in copy number data
cn1 = read.table("data/copy_number_tum1")
cn1 = cn1[,c(1,2,3,5)]
cn2 = read.table("data/copy_number_tum2")
cn2 = cn2[,c(1,2,3,5)]
cn3 = read.table("data/copy_number_tum3")
cn3 = cn3[,c(1,2,3,5)]

#set sample names
names = c("Sample1","Sample2","Sample3")


#make an output directory, deleting old results first if they exist
suppressWarnings(dir.create("results"))
unlink("results/*", recursive=TRUE)


print("")
print("=========================================================")
print("Test 1 - single sample - shortTest1")
print("")
#run one sample
sciClone(vafs=v1,
         copyNumberCalls=cn1,
         sampleNames=names[1],
         outputPrefix="results/shortTest1",
         overlayClusters=TRUE)




## #run only one sample, but all sites are removed by excluded regions
## #should fail with "can't do clustering - no copy number 2 regions to operate on in sample 1"
## sciClone(vafs=v1,
##          regionsToExclude=regions,
##          copyNumberCalls=cn1,
##          sampleNames=names,
##          outputPrefix="test.results/shortTest1",
##          overlayClusters=TRUE)

print("")
print("=========================================================")
print("Test 2 - two samples - shortTest2")
print("")
#run two samples
sciClone(vafs=list(v1,v2),
         copyNumberCalls=list(cn1,cn2),
         sampleNames=names[1:2],
         outputPrefix="results/shortTest2",
         overlayClusters=TRUE)



print("")
print("=========================================================")
print("Test 2 - three samples - shortTest3")
print("")
#run two samples
sciClone(vafs=list(v1,v2,v3),
         copyNumberCalls=list(cn1,cn2,cn3),
         sampleNames=names,
         outputPrefix="results/shortTest3",
         overlayClusters=TRUE)


