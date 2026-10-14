use strict;
use warnings;

if($#ARGV != 3)
{
  print "ERROR: There should be 4 parameters.\n";
  print "perl getPdbqt.pl 0 1 2 3\n";
  print "0: pdbqt output file from docking\n";
  print "1: output log file from docking\n";
  print "2: prefix for pdbqt files which will be generated\n";
  print "3: frame skip\n";
  exit 1;
}

my $fileIn1    = $ARGV[0];
my $fileIn2    = $ARGV[1];
my $fileOutPre = $ARGV[2];
my $skip       = $ARGV[3];
my $fileOut;
my $line;
my @words;
my @file1data;
my @file2data1;
my @file2data2;
my @file2data3;
my $dataNum1;
my $dataNum2;
my $dataNum3;
my $dataCount2;
my $dataCount3;
my $scale = 0.375000;

open(FIN, "$fileIn1") or die "ERROR: Cannot open file $fileIn1\n";
while($line = <FIN>)
{
  chomp($line);
  push(@file1data, "$line\n");
}
close FIN;

open(FIN, "$fileIn2") or die "ERROR: Cannot open file $fileIn2\n";
$dataNum1 = 0;
$dataNum2 = 0;
$dataNum3 = 0;
$dataCount2 = 0;
$dataCount3 = 0;
while($line = <FIN>)
{
  if($line =~ /^COOR/)
  {
    my $lineNum1 = 0;
    my $lineNum2 = 0;
    my $lineNum3 = 0;
    while($line = <FIN>)
    {
      if($line =~ /^ATOM1/)
      {
        @words = split(/\s+/, $line);
        push(@{$file2data1[$dataNum1][$lineNum1]}, $words[2] * $scale);
        push(@{$file2data1[$dataNum1][$lineNum1]}, $words[3] * $scale);
        push(@{$file2data1[$dataNum1][$lineNum1]}, $words[4] * $scale);
#         print("$lineNum\n");
        $lineNum1++;
      }
      elsif($line =~ /^ATOM2/)
      {
        if($dataCount2 == 0)
        {
          @words = split(/\s+/, $line);
          push(@{$file2data2[$dataNum2][$lineNum2]}, $words[2] * $scale);
          push(@{$file2data2[$dataNum2][$lineNum2]}, $words[3] * $scale);
          push(@{$file2data2[$dataNum2][$lineNum2]}, $words[4] * $scale);
#           print("$lineNum\n");
          $lineNum2++;
        }
      }
      elsif($line =~ /^ATOM3/)
      {
        if($dataCount3 == 0)
        {
          @words = split(/\s+/, $line);
          push(@{$file2data3[$dataNum3][$lineNum3]}, $words[2] * $scale);
          push(@{$file2data3[$dataNum3][$lineNum3]}, $words[3] * $scale);
          push(@{$file2data3[$dataNum3][$lineNum3]}, $words[4] * $scale);
#           print("$lineNum\n");
          $lineNum3++;
        }
      }
      elsif($line =~ /^DONE1/)
      {
        print("$line\n");
        $dataNum1++;
        last;
      }
      elsif($line =~ /^DONE2/)
      {
        print("$line\n");
        $dataCount2++;
        if(($dataCount2 % $skip) == 0)
        {
          $dataNum2++;
          $dataCount2 = 0;
        }
        last;
      }
      elsif($line =~ /^DONE3/)
      {
        print("$line\n");
        $dataCount3++;
        if(($dataCount3 % $skip) == 0)
        {
          $dataNum3++;
          $dataCount3 = 0;
        }
        last;
      }
    }
  }
}
close FIN;

for(my $i = 0; $i < $dataNum1; $i++)
{
  my $lineNum = $#{$file2data1[$i]} + 1;
  for(my $j = 0; $j < $lineNum; $j++)
  {
    print("$i $j $file2data1[$i][$j][0] $file2data1[$i][$j][1] $file2data1[$i][$j][2]\n");
  }

  $fileOut = "${fileOutPre}_1_$i.pdbqt";
  open(FOUT, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  $fileOut = "${fileOutPre}_1_$i.pdb";
  open(FOUT2, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  my $atomIndex = 0;
  foreach my $line (@file1data)
  {
    if($line =~ /^ATOM/)
    {
      my $line2 = $line;
      substr($line2, 30, 8, sprintf("%8.3f", $file2data1[$i][$atomIndex][0]));
      substr($line2, 38, 8, sprintf("%8.3f", $file2data1[$i][$atomIndex][1]));
      substr($line2, 46, 8, sprintf("%8.3f", $file2data1[$i][$atomIndex][2]));
      print FOUT "$line2";
      print FOUT2 "$line2";
      $atomIndex++;
    }
    else
    {
      print FOUT "$line";
    }
  }
  close FOUT;
  close FOUT2;
}

for(my $i = 0; $i < $dataNum2; $i++)
{
  my $lineNum = $#{$file2data2[$i]} + 1;
  for(my $j = 0; $j < $lineNum; $j++)
  {
    print("$i $j $file2data2[$i][$j][0] $file2data2[$i][$j][1] $file2data2[$i][$j][2]\n");
  }

  $fileOut = "${fileOutPre}_2_$i.pdbqt";
  open(FOUT, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  $fileOut = "${fileOutPre}_2_$i.pdb";
  open(FOUT2, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  my $atomIndex = 0;
  foreach my $line (@file1data)
  {
    if($line =~ /^ATOM/)
    {
      my $line2 = $line;
      substr($line2, 30, 8, sprintf("%8.3f", $file2data2[$i][$atomIndex][0]));
      substr($line2, 38, 8, sprintf("%8.3f", $file2data2[$i][$atomIndex][1]));
      substr($line2, 46, 8, sprintf("%8.3f", $file2data2[$i][$atomIndex][2]));
      print FOUT "$line2";
      print FOUT2 "$line2";
      $atomIndex++;
    }
    else
    {
      print FOUT "$line";
    }
  }
  close FOUT;
  close FOUT2;
}

for(my $i = 0; $i < $dataNum3; $i++)
{
  my $lineNum = $#{$file2data3[$i]} + 1;
  for(my $j = 0; $j < $lineNum; $j++)
  {
    print("$i $j $file2data3[$i][$j][0] $file2data3[$i][$j][1] $file2data3[$i][$j][2]\n");
  }

  $fileOut = "${fileOutPre}_3_$i.pdbqt";
  open(FOUT, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  $fileOut = "${fileOutPre}_3_$i.pdb";
  open(FOUT2, "> $fileOut") or die "ERROR: Cannot open file $fileOut\n";
  my $atomIndex = 0;
  foreach my $line (@file1data)
  {
    if($line =~ /^ATOM/)
    {
      my $line2 = $line;
      substr($line2, 30, 8, sprintf("%8.3f", $file2data3[$i][$atomIndex][0]));
      substr($line2, 38, 8, sprintf("%8.3f", $file2data3[$i][$atomIndex][1]));
      substr($line2, 46, 8, sprintf("%8.3f", $file2data3[$i][$atomIndex][2]));
      print FOUT "$line2";
      print FOUT2 "$line2";
      $atomIndex++;
    }
    else
    {
      print FOUT "$line";
    }
  }
  close FOUT;
  close FOUT2;
}

exit;
