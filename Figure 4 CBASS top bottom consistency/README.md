Simple R script for analysis of CBASS ED50 top/bottom assignment consistency;
user can chosse the number of top/bottom colonies to consider assigning X in "item_per_group <- X" (line 22);
reads data and checks consistency/correspondence of top/bottom assignment based on ED50s of two (consecutive) CBASS runs; 
creates subsets of data to report on consistency in dependency of number of colonies assayed and numbger of top/bottom assignments considered (e.g. top/bottom5 of 20 assessed colonies)
infile is csv with ID, ED50_1, ED50_2; additiional columnes will be ignored
