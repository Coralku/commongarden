Simple R script for analysis of CBASS ED50 top/bottom assignment consistency;

User can chosse the number of top/bottom colonies to consider by assigning X in "item_per_group <- X" (line 22);

Sript reads data and checks consistency/correspondence of top/bottom assignment based on ED50s of two (consecutive) CBASS runs; 

Script creates subsets of data to report on consistency in dependency of number of colonies assayed and numbger of top/bottom assignments considered (e.g. top/bottom5 of 20 assessed colonies);

Infile is csv with ID, ED50_1, ED50_2; additional columns will be ignored
