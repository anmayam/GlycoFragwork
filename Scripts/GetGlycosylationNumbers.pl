# Anoop Mayampurath
# Program to parse a glycopeptide list and 
# get numbers such as number of glycoproteins, glycosylational sites
# and glycopeptides

#use strict ; 
#use warnings ;
#
# MAKE SURE tHAT INPUT IS SORTED ACCORDING TO GLYCAN, Site AND PROTEIN 
#
# # V3 version looks at site instead of peptide
#  Format of input should be
#  ID, Protein, Sequence, Site, Glycan_composition (includes headers)

my $glycoproteinList = $ARGV[0] ; 
my @temp ; 
my @temp1 ; 
my @temp2 ; 
my $line ; 
open(IN, "<$glycoproteinList") or die ; 

my %gp_list ; 
my @gp_list_values; 
my %gp_num_sites ; 
my %gp_num_glycopeptides ; 
my %gp_unique_sites ; 

my $total_gp = 0 ; 
my $total_sites = 0 ; 
my $total_glycopeptides = 0 ; 
my $i ; 
my $j ; 
my $found_site ; 
my $found_glycopeptide ; 

my $peptide ; 
my $intactpeptide ; 
my $glycopeptide ; 
my $glycan ; 
my $this_site ; 
my $site ; 

# Read off header
$line = <IN> ; 
my $break = 1; 

# Start
while ($line = <IN>)
{
	chomp($line) ; 
	@temp = split(',', $line ) ;
        if ($temp[0] eq 4)
	{
		$break = 1; 
	}
	if ($temp[0] eq 13273)
	{
		$break = 0 ; 
	}	
	if (exists $gp_list{$temp[1]})
	{
		# Gp exists, see if the site has been picked up
		for (keys %gp_list)
		{
			if ($_ eq $temp[1])
			{				
				# Get all the glycopeptides out so far
				@gp_list_values = @{$gp_list{$_}} ; 
				$found_site = 0 ; 
				$found_glycopeptide = 0 ; 
				for ($i = 0 ; $i < scalar(@gp_list_values); $i++)
				{
					$glycopeptide = $gp_list_values[$i] ;
					@temp1 = split(',', $glycopeptide) ; 				
					$site = $temp1[0] ; 
					$glycan = $temp1[1] ; 

					if (($site =~ m/\//) or ($temp[3] =~ m/\//))
					{

					}
					if ($temp[3] =~ m/$site/)
					{
						$found_site = 1 ; 
						if ($glycan eq $temp[4])
						{		
							#glycan is already present 
							$found_glycopeptide = 1 ; 
						}						
					}
				}								
			        if (!$found_site)
				{
					#	print "$temp[0], $temp[2]\n" ; 
					$total_sites = $total_sites +1 ;
				        $this_site = $temp[1].",".$temp[3] ; 
					$gp_unique_sites{$temp[0]} = $this_site ; # Add site 
				}
				if (!$found_glycopeptide)
				{
					# same site, new glycan.
					$total_glycopeptides = $total_glycopeptides +1 ; 
					$glycopeptide = $temp[3].",".$temp[4] ; 
					push(@gp_list_values, $glycopeptide) ; 
					#$gp_list{$_} = \@gp_list_values ; 
					push (@{ $gp_list{$_}}, $glycopeptide) ; 
				}
			   }
		   }
	   } 
	   else
	   {
		   # new glycoprotein
		   $total_gp = $total_gp +1 ; 
		   $total_glycopeptides = $total_glycopeptides + 1; 
		   # print "$temp[0],$temp[2]\n" ;
		   $total_sites = $total_sites + 1; 


		   @gp_list_values = () ; 
		   $glycopeptide = $temp[3].",".$temp[4] ; 
		   $gp_list{$temp[1]} = [$glycopeptide];

		   # add site
	   	   $this_site = $temp[1].",".$temp[3] ; 
		   $gp_unique_sites{$temp[0]} = $this_site ;   

	   }
			
}


print "Number of glycoproteins : $total_gp\n" ; 
print "Number of glycosylation sites : $total_sites\n" ; 
print "Number of glycopeptides: $total_glycopeptides\n" ; 

