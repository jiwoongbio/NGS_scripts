sub getCombinationList {
	my (@tokenListList) = @_;
	my @combinationList = ();
	if(my ($index) = grep {ref($tokenListList[$_])} 0 .. $#tokenListList) {
		foreach(@{$tokenListList[$index]}) {
			push(@combinationList, getCombinationList(@tokenListList[0 .. ($index - 1)], $_, @tokenListList[($index + 1) .. $#tokenListList]));
		}
	} else {
		push(@combinationList, \@tokenListList);
	}
	return @combinationList;
}

sub getBitNumberList {
	my ($number) = @_;
	my @bitList = ();
	while($number) {
		push(@bitList, $number & 1);
		$number >>= 1;
	}
	return map {2 ** $_} grep {$bitList[$_]} 0 .. $#bitList;
}
