use strict;
use warnings;
# v.0.0.1

if($#ARGV != 2)
{
  print("ERROR: The number of required argument is 3, but there are ", $#ARGV + 1, " argument(s) read.\n\n");

  print("USAGE: perl $0 1 2 3\n");
  print("  1: the input pdbqt file\n");
  print("  2: the output pdbqt file\n");
  print("  3: the file with atom pairs information for distance fixing\n");
  print("\n");
  exit;
}

my $pdbqtIn    = "$ARGV[0]";
my $pdbqtOut   = "$ARGV[1]";
my $atomPairIn = "$ARGV[2]";
my $fileIn;
my $fileOut;
my @pdbqtData;
my $pdbqtDataLen;
my @atomMarker;
my @atomPairList;
my $atomPairListLen;
my $line;
my @words;
my $oriAtomNum;
# my @oriAtomSerial;
my @newAtomSerial;
my $tmpStr;

$fileIn = "$ARGV[0]";
open(FIN, "$fileIn") or die "ERROR: Cannot open file $fileIn.\n";
while($line = <FIN>)
{
  chomp($line);
  push(@pdbqtData, $line);
  if($line =~ /^ATOM|^HETATM/)
  {
    push(@newAtomSerial, int(substr($line, 6, 5)));
    $tmpStr  = substr($line, 12, 4);
    $tmpStr .= " " . substr($line, 17, 3);
    $tmpStr .= " " . substr($line, 21, 1);
    $tmpStr .= " " . substr($line, 22, 4);
    $tmpStr .= " " . substr($line, 77, 1);
    $tmpStr =~ s/^\s+|\s+$//g;
    $tmpStr =~ s/\s+/ /g;
    push(@atomMarker, $tmpStr);
  }
}
close FIN;
$oriAtomNum = $#atomMarker + 1;
print("There are $oriAtomNum atoms in file $pdbqtIn.\n");

$fileIn = "$ARGV[2]";
open(FIN, "$fileIn") or die "ERROR: Cannot open file $fileIn.\n";
while($line = <FIN>)
{
  chomp($line);
  @words = split(/\:/, $line);
  $words[0] =~ s/^\s+|\s+$//g; # remove leading and trailing white spaces
  $words[1] =~ s/^\s+|\s+$//g;
  $words[0] =~ s/\s+/ /g; # convert multiple consecutive spaces into one space
  $words[1] =~ s/\s+/ /g;
  push(@{$atomPairList[0]}, $words[0]);
  push(@{$atomPairList[1]}, $words[1]);
}
close FIN;

$atomPairListLen = $#{$atomPairList[0]} + 1;
$pdbqtDataLen  = $#pdbqtData + 1;

print("$atomPairListLen atom pair(s) read:\n");
for(my $atomPairListIndex = 0; $atomPairListIndex < $atomPairListLen; $atomPairListIndex++)
{
  print("  (", $atomPairListIndex + 1, ") $atomPairList[0][$atomPairListIndex] : $atomPairList[1][$atomPairListIndex]\n");
}

# Start to get the mapping from the original serials to the new serials
for(my $atomPairListIndex = 0; $atomPairListIndex < $atomPairListLen; $atomPairListIndex++)
{
  my $atomIndex = -1;
  for(my $dataLineIndex = 0; $dataLineIndex < $pdbqtDataLen; $dataLineIndex++)
  {
    if($pdbqtData[$dataLineIndex] =~ /^ATOM|^HETATM/)
    {
      $atomIndex++;
      my $foundPair = 0;
      for(my $i = 0; $i < 2; $i++)
      {
        if($atomPairList[$i][$atomPairListIndex] eq $atomMarker[$atomIndex])
        {
          $foundPair++;
          my $i2 = ($i + 1) % 2;

          for(my $aIndex = $atomIndex + 1; $aIndex < $oriAtomNum; $aIndex++)
          {
            $newAtomSerial[$aIndex]++;
          }

          my $atomIndex2 = $atomIndex;
          for(my $dataLineIndex2 = $dataLineIndex + 1; $dataLineIndex2 < $pdbqtDataLen; $dataLineIndex2++)
          {
            if($pdbqtData[$dataLineIndex2] =~ /^ATOM|^HETATM/)
            {
              $atomIndex2++;
              if($atomPairList[$i2][$atomPairListIndex] eq $atomMarker[$atomIndex2])
              {
                $foundPair++;

                for(my $aIndex = $atomIndex2 + 1; $aIndex < $oriAtomNum; $aIndex++)
                {
                  $newAtomSerial[$aIndex]++;
                }

                last;
              }
            }
          }

          last;
        }
      }
      if($foundPair == 1)
      {
        print("ERROR: Cannot find the same atom pair \"$atomPairList[0][$atomPairListIndex] : $atomPairList[1][$atomPairListIndex]\" in file $pdbqtIn.\n");
        exit;
      }
      if($foundPair == 2)
      {
        last;
      }
    }
  }
}
# End of getting the mapping the original serials to the new serials


# Start to add and modify atoms
for(my $atomPairListIndex = 0; $atomPairListIndex < $atomPairListLen; $atomPairListIndex++)
{
  print("Add and modify atoms for atom pair \"(", $atomPairListIndex + 1, ") $atomPairList[0][$atomPairListIndex] : $atomPairList[1][$atomPairListIndex]\"\n");
  my $atomIndex = -1;
  for(my $dataLineIndex = 0; $dataLineIndex < $pdbqtDataLen; $dataLineIndex++)
  {
    if($pdbqtData[$dataLineIndex] =~ /^ATOM|^HETATM/)
    {
      $atomIndex++;
      my $foundPair = 0;
      for(my $i = 0; $i < 2; $i++)
      {
        if($atomPairList[$i][$atomPairListIndex] eq $atomMarker[$atomIndex])
        {
          $foundPair++;
          my $i2 = ($i + 1) % 2;
          $dataLineIndex++;
          splice(@pdbqtData, $dataLineIndex, 0, "ATOM       ");
          $pdbqtDataLen++;
          $atomIndex++;
          splice(@atomMarker, $atomIndex, 0, "");

          my $atomIndex2 = $atomIndex;
          for(my $dataLineIndex2 = $dataLineIndex + 1; $dataLineIndex2 < $pdbqtDataLen; $dataLineIndex2++)
          {
            if($pdbqtData[$dataLineIndex2] =~ /^ATOM|^HETATM/)
            {
              $atomIndex2++;
              if($atomPairList[$i2][$atomPairListIndex] eq $atomMarker[$atomIndex2])
              {
                $foundPair++;
                splice(@pdbqtData, $dataLineIndex2 + 1, 0, "ATOM       ");
                $pdbqtDataLen++;
                $atomIndex2++;
                splice(@atomMarker, $atomIndex2, 0, "");

                substr($pdbqtData[$dataLineIndex-1], 77, 6, substr($atomPairList[$i][$atomPairListIndex], -1, 1) . sprintf("G%-4s", $atomPairListIndex));
                $pdbqtData[$dataLineIndex] = $pdbqtData[$dataLineIndex2];
                my $aIndex = $atomIndex;
                for(my $dataLineIndex3 = $dataLineIndex; $dataLineIndex3 < $pdbqtDataLen; $dataLineIndex3++)
                {
                  if($pdbqtData[$dataLineIndex3] =~ /^ATOM|^HETATM/)
                  {
                    my $atomSerial = $aIndex + 1;
                    substr($pdbqtData[$dataLineIndex3], 6, 5, sprintf("%5s", $atomSerial));
                    $aIndex++;
                  }
                }
                @words = split(/\s/, $atomPairList[$i2][$atomPairListIndex]);
                substr($pdbqtData[$dataLineIndex], 12, 4, sprintf("%4s", $words[4] . "GX"));
                @words = split(/\s/, $atomPairList[$i][$atomPairListIndex]);
                substr($pdbqtData[$dataLineIndex], 17, 3, sprintf("%3s", $words[1]));
                substr($pdbqtData[$dataLineIndex], 22, 4, sprintf("%4s", $words[3]));
                substr($pdbqtData[$dataLineIndex], 77, 6, sprintf("G%-5s", $atomPairListIndex));

                substr($pdbqtData[$dataLineIndex2], 77, 6, substr($atomPairList[$i2][$atomPairListIndex], -1, 1) . sprintf("G%-4s", $atomPairListIndex));
                $dataLineIndex--;
                $dataLineIndex2++;
                $pdbqtData[$dataLineIndex2] = $pdbqtData[$dataLineIndex];
                $aIndex = $atomIndex2;
                for(my $dataLineIndex3 = $dataLineIndex2; $dataLineIndex3 < $pdbqtDataLen; $dataLineIndex3++)
                {
                  if($pdbqtData[$dataLineIndex3] =~ /^ATOM|^HETATM/)
                  {
                    my $atomSerial = $aIndex + 1;
                    substr($pdbqtData[$dataLineIndex3], 6, 5, sprintf("%5s", $atomSerial));
                    $aIndex++;
                  }
                }
                @words = split(/\s/, $atomPairList[$i][$atomPairListIndex]);
                substr($pdbqtData[$dataLineIndex2], 12, 4, sprintf("%4s", $words[4] . "GX"));
                @words = split(/\s/, $atomPairList[$i2][$atomPairListIndex]);
                substr($pdbqtData[$dataLineIndex2], 17, 3, sprintf("%3s", $words[1]));
                substr($pdbqtData[$dataLineIndex2], 22, 4, sprintf("%4s", $words[3]));
                substr($pdbqtData[$dataLineIndex2], 77, 6, sprintf("G%-5s", $atomPairListIndex));
                last;
              }
            }
          }

          last;
        }
      }
      if($foundPair == 1)
      {
        print("ERROR: Bug?! Cannot find the same atom pair \"$atomPairList[0][$atomPairListIndex] : $atomPairList[1][$atomPairListIndex]\" in file $pdbqtIn.\n");
        exit;
      }
      if($foundPair == 2)
      {
        last;
      }
    }
  }
}
# End of adding and modifying atoms


# Appyly new atom serial to rows of BRANCH and ENDBRANCH
for(my $dataLineIndex = 0; $dataLineIndex < $pdbqtDataLen; $dataLineIndex++)
{
  if($pdbqtData[$dataLineIndex] =~ /^BRANCH|^ENDBRANCH/)
  {
    @words = split(/\s+/, $pdbqtData[$dataLineIndex]);
    $words[1] = $newAtomSerial[$words[1] - 1];
    $words[2] = $newAtomSerial[$words[2] - 1];
    $pdbqtData[$dataLineIndex] = join(" ", @words);
  }
}

$fileOut = "$ARGV[1]";
open(FOUT, "> $fileOut") or die "ERROR: Cannot open file $fileOut.\n";
for(my $dataLineIndex = 0; $dataLineIndex < $pdbqtDataLen; $dataLineIndex++)
{
  print FOUT "$pdbqtData[$dataLineIndex]\n";
}
close FIN;

exit;
