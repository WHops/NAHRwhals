# T2T is denoted with cryptic chromsome names. I'm sure there are good reasons
# to do that, but for readability, I prefer to display simple chromosome names
# in output plots involving T2T. This function helps with that. Give it a cryptic
# name and it returns a simple one. If cryptic name is unknown, return the original
# name. String in, string out. 
#' @export
translate_t2t_chr_to_readable <- function(chrname){
  
  t2t_chrs_translate = data.frame(
    'CP068277.2' =	'chr1',
    'CP068276.2' =	'chr2',
    'CP068275.2' =	'chr3',
    'CP068274.2' =	'chr4',
    'CP068273.2' =	'chr5',
    'CP068272.2' =	'chr6',
    'CP068271.2' =	'chr7',
    'CP068270.2' =	'chr8',
    'CP068269.2' =	'chr9',
    'CP068268.2' =	'chr10',
    'CP068267.2' =	'chr11',
    'CP068266.2' =	'chr12',
    'CP068265.2' =	'chr13',
    'CP068264.2' =	'chr14',
    'CP068263.2' =	'chr15',
    'CP068262.2' =	'chr16',
    'CP068261.2' =	'chr17',
    'CP068260.2' =	'chr18',
    'CP068259.2' =	'chr19',
    'CP068258.2' =	'chr20',
    'CP068257.2' =	'chr21',
    'CP068256.2' =	'chr22',
    'CP068255.2' =	'chrX',
    'CP086569.2' =	'chrY'
  )

  # Return the translation, or orig if translation is impossible: 
  if (chrname %in% colnames(t2t_chrs_translate)){
    return(as.character(t2t_chrs_translate[chrname]))
  } else { # If input name is not a contig name in T2T (ideally this is never called...)
    print("Chromsosome name unknown. Not attemptying to translate name.")
    return(chrname)
  }
  
  
}
