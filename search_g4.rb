#This script is used to search for G4 motifs across Parhyale hawaiensis genome using the G4Hunter algorithm
#It requires to run on conda py2 enviorenment

#This script is used to search for G4 motifs across Parhyale hawaiensis genome using the G4Hunter algorithm
puts "Scaffold" + "\t" + "Start" + "\t" +  "End" + "\t" + "Sequence" + "\t" + "Length" + "\t" + "Score" + "\t" + "NBR"
file_path= "/drives/ssd1/manuel/phaw/2021_analysis/annotation/utr_analysis/phaw_5_utrs.fasta.tbl"
file_name= file_path.split("/")[-1].split(".")[0]
#File.open("parhyale_hawaiensis_15Feb20182_1IWOV.fa.tbl","r").each do |line| #Whole-genome search

File.open(file_path,"r").each do |line|
    
    scaffold = []
    #File.open("tempfile.fasta","w") do |f|
    File.open("#{file_name}.fasta","w") do |f|
        id= line.split("\t")[0] 
        seq= line.chomp.split("\t")[1]
        f.puts ">" + id
        f.puts seq
        scaffold << id
    end
 `python2.7 /drives/raid/AboobakerLab/software/G4-hunter/G4Hunter.py -i "#{file_name}.fasta" -o . -w 25 -s 0.9`

     #File.open("./Results_tempfile/tempfile-Merged.txt","r").each_line do |o|
     File.open("./Results_#{file_name}/#{file_name}-Merged.txt","r").each_line do |o|
         #puts o
        if o[0] =~  /[[:digit:]]/
            #if (o.split("\t")[4].to_f).abs >= 1 #Print G4s found in both strands
            if o.split("\t")[4].to_f >= 1 #Print G4s found in sense strands
			puts scaffold[0] + "\t" + o
            end
         end
     end
  
end
#`rm tempfile.fasta`
`rm "#{file_name}.fasta"`
