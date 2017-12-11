import os
import periodictable as pt

tagdir= os.getenv('TAGSPACE_DATA')

numbers = [elem.number for elem in pt.elements]
symbols = [elem.symbol for elem in pt.elements]

ptnumdict = dict(zip(symbols,numbers))
ptsymdict = dict(zip(numbers,symbols))