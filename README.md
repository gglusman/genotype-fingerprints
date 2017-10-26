# genotype-fingerprints
Software for creating and comparing genotype fingerprints.  
More information and datasets: http://db.systemsbiology.net/gestalt/genotype_fingerprints/  
Preprint: TBD

1. To create a fingerprint for a genome:  
	`bin/computeGF.pl` _myGenome path-to-my-genotypes/myGenome.gz_  
	...will generate myGenome.outn.gz (and some other files)

2. To compare two fingerprints:  
	`bin/compareGFs.pl` _myFirstGenome.outn.gz mySecondGenome.outn.gz_

3. To serialize fingerprints into a database, using L=1000:  
	`bin/serializeGFs.pl` _myFingerprintCollection_ 1000 _@myListOfFingerprints_  
	`bin/serializeGFs.pl` _myFingerprintCollection_ 1000 _*.outn.gz_

4. To compare a fingerprint to a database:  
	`bin/searchGFs.pl` _myFirstGenome.outn.gz myFingerprintCollection_  
	...see the data directory for an example database (Corpas family)

5. To compare two databases:  
	`bin/searchGFs.pl` _aFingerprintCollection anotherFingerprintCollection_

6. To perform all-against-all comparisons in one database:  
	`bin/searchGFs.pl` _aFingerprintCollection_

