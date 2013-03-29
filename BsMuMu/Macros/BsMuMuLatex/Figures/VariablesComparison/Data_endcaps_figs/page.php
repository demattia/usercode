<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html>

<head>
  <title>Data Endcaps</title>
</head>

<body>
  <h3>Data Endcaps</h3>
  <?
     $handle = opendir(dirname(realpath(__FILE__)).'');
     while($file = readdir($handle)) {
     if($file !== '.' && $file !== '..') {
     echo '<img src="'.$file.'" border="0" />';
     }
     }
  ?>
</body>

</html>

