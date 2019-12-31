sub doSys #( string command )
{
    my $systemCommand = $_[0];
    
    system("date");
    print STDERR "$0: Start to execute   [$systemCommand] \n";

    #my $returnCode = system( "echo test" );    
    my $returnCode = system( $systemCommand );
    
    if ( $returnCode != 0 ) 
    { 
        print STDERR "$0: Failed executing     [$systemCommand] with doSys exit code $returnCode\n";
        system("date");
        exit ($returnCode);  
    }
    else
    {
        print STDERR "$0: Finished executing [$systemCommand] with doSys exit code $returnCode\n";
        system("date");
        return 0;
    }
    
}

sub doSys2
{
	my $systemCommand = $_[0];
	system("date");
	print STDERR "$0: Start to execute   [$systemCommand] \n";
	my $returnCode = system( $systemCommand )>>8; ## perl exit code has been multiplied by 256, thus revert back
	print STDERR "$0: Finished executing [$systemCommand] with doSys exit code $returnCode\n";
	system("date");
	return $returnCode;
}

1;
