<html>
 <head>
  <title>Journal</title>
 </head>
 <body>

<h1> Personal journal of Marco De Mattia </h1>
<div>Disclaimer: If I did not send you the link to a specific page you can assume the information on the links is typically work in progress.</div>

<?php

$path = realpath('./');
$directories;
$date;

$objects = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($path), RecursiveIteratorIterator::SELF_FIRST);
foreach($objects as $name){
  if( preg_match("/html$/", $name) ) {
    $linkString = substr($name, 40);
    $splitted = explode("/", $linkString);
    $i = 0;
    $numDirs = count($splitted) - 1;
    while( $numDirs > $i ) {
      if( $directories[$i] != $splitted[$i] ) {
	$directories[$i] = $splitted[$i];
	if($i == 1) {
	  echo "<h2> $splitted[$i] </h3>";
	}
	if($i == ($numDirs - 1)) {
	  $date = $splitted[$i];
	}
      }
      $i++;
    }
    # echo "<a href=./$linkString>$linkString</a><br>";
    echo "<a href=./$linkString>$date</a><br>";
  }
}

?>

 </body>
</html>
