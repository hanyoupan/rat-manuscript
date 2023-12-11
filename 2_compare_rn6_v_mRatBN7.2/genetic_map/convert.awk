{
	s = $2 "\t" $3
	for (i=6;i<=NF;i+=3){
		max=$i
		if ($(i+1) > max)
			max = $(i+1)
                if ($(i+2) > max)
                        max = $(i+2)
		s = s "\t" ($i / max) " " ($(i+1) / max) " 0 0 " ($(i+2) / max) " 0 0 0 0 0"
	}
	print s
}
