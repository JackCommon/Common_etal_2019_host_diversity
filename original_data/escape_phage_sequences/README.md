These are all the sequence (.seq) and chromatogram (.ab1) files returned from
Sanger sequencing (SourceBiosciences, Sheffield, UK). I haven't yet gotten 
around to giving them useful names, so you can read the files as follows:

	- The numbers and "jc" at the start can be safely ignored
	- The next section is which phage it is:
		- "Ancestral" = ancestral DMS3vir
		- "Phage" = pre-evolved phage from the library
		- The timepoint (Tn), replicate (Rn), phage isolate (1-12)
	- The sequencing primer used, which can be found in the Supp Info

Hope that all makes sense and one day I'll get round to writing a regex to 
batch rename these awkward files! But to quote Aragorn, "it is not this day!"
