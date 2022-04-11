reset session
set term pngcairo

FILES = system("ls ../data/*.txt")      # Linux
#FILES = system("dir /B *.dat")  # Windows

myOutput(s) = sprintf("%s_3f.png",s)

do for [FILE in FILES] {
    set output myOutput(FILE)
    plot FILE u 1:9 w lines lc rgb "navy" title FILE
}
set output