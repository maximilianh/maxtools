fastaDistribute:
	rm -rf outDir
	mkdir outDir
	fastaDistribute fastaDistribute.fa outDir
	cat outDir/*
	ls -la outDir | grep -v ^d | cut -d' ' -f7,10 > fastaDistribute.out
	cmp fastaDistribute.out fastaDistribute.expect
