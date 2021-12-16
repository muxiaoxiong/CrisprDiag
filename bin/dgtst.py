import primer3
def test():
	lst=["A","T","C","G"]
	n=0
	word=""
	for a in lst:
		for b in lst:
			for c in lst:
				for d in lst:
					for e in lst:
						for f in lst:
							word+=a+b+c+d+e+f
							if dg(word) <= -4:
								n+=1
								print(dg(word))
								print(word)
							word=""
	print(n)
def dg(seq):
    flag=seq[-6:]
    dg=primer3.calcHomodimer(flag).dg
    print(dg)
    return dg

dg("TAGCTA")