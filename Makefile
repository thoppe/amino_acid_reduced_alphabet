bead_targets = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
bead_targets = 1 2
cmd = python enumerate_alphabets.py

all:
	echo "pass"

MJ96:
	$(foreach b, $(bead_targets), $(cmd) MJ96 -b $(b);)

BT:
	$(foreach b, $(bead_targets), $(cmd) BT -b $(b);)

