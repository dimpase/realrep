Read("interesting_cases.g");
SetAssertionLevel(10);
g:=GL(3,2);
t:=CharacterTable(g);
irr:=Irr(t);
m:=make_real_representation_NC(g, irr[4]);
Print(List(m.generators, x-> x=ComplexConjugate(x)), "\n should be [true,true]\n");
