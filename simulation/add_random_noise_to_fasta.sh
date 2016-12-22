# add_random_noise_to_fasta.sh - Deos exactly what it says on the tin!
# Assumes that the sequences have been generated from seq-gen and transformed to fasta using readseq, which leaves the names as ">TAXANAME LENGTH bp"
# if the statement "int(rand()*100)" then randomly 1% of sites will be changed randomly, to change this to 10% change this to int(rand()*10)
# 18 Aug 2016
# One Argument needed, the FASTA file. In the script, $1 is the fasta file, suvh as "test_tree1_L1000.nex.fasta"
# This will randomly select an amino acid to change the base to, it doesn not follow any models (for now)

if [ "$#" -ne 1 ]; then
head -6 add_random_noise_to_fasta.sh
else

random=$RANDOM;
sed '/>/s/$/	/g' < $1 | tr -d '\n' | sed 's/>/\
>/g'  | awk -v RAND=$random 'BEGIN{srand(RAND); aa[0] = "A";aa[1] = "C";aa[2] = "D";aa[3] = "E";aa[4] = "F";aa[5] = "G";aa[6] = "H";aa[7] = "I";aa[8] = "K";aa[9] = "L";aa[10] = "M";aa[11] = "N";aa[12] = "P";aa[13] = "Q";aa[14] = "R";aa[15] = "S";aa[16] = "T";aa[17] = "V";aa[18] = "W";aa[19] = "Y";;};{ print $1; split($4, string, ""); for(i=1; i<=length($4); i++) { if (int(rand()*100) == 0 ) {printf aa[int(rand()*20)];} else {printf string[i]}}; printf "\n";}'

fi