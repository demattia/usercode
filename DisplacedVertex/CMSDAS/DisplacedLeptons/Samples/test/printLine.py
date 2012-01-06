file = open("runAllAnalysis.sh")
out = open("out.sh", "w")
for line in file:
    out.write("echo \""+line.rstrip()+"\"\n")
    out.write(line)
