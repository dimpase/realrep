Read("interesting_cases.g");
SetAssertionLevel(10);
g:=GL(3,2);
t:=CharacterTable(g);
irr:=Irr(t);
m:=make_real_representation_NC(g, irr[4]);
