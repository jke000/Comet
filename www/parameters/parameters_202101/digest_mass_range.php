<?php include "head.php" ; ?>
<html>
<body>
<?php include "topmenu.php" ; ?>
<?php include "imgbar.php" ; ?>

<div id="page">
   <div id="content_full">
      <div class="post hr">

         <h2>Comet parameter: digest_mass_range</h2>

         <ul>
         <li>Defines the mass range of peptides to search (based on MH+ or the singly
         protonated mass).
         <li>This parameter has two decimal values.
         <li>The first value is the lower mass cutoff and the second value is
         the high mass cutoff.
         <li>Only spectra with experimental MH+ masses within (or equal to) the defined
         mass ranges are searched.
         <li>Valid values are two decimal numbers where the first number must
         be less or equal to the second number.
         <li>The default value is "600.0 8000.0" if this parameter is missing.
         </ul>

         <p>Example:
         <br><tt>digest_mass_range = 600.0 8000.0</tt> &nbsp; &nbsp; <i>search only 600.0 to 8000.0 mass range</i>
         <br><tt>digest_mass_range = 400.0 5000.0</tt> &nbsp; &nbsp; <i>search all spectra with peptide masses between 400.0 and 5000.0</i>

      </div>
   </div>
   <div style="clear: both;">&nbsp;</div>
</div>

<?php include "footer.php" ; ?>
</html>