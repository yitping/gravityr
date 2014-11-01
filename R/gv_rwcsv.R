#' Read and write a custom formatted CSV file
#' 
#' Such a file stores multiple 1D/2D arrays in a csv format. Each array is
#' preceeded with a header describing its name, row and colum. This format is
#' currently used for storing gvspc calibration data.
#' 
#' @param filein Path and name of the csv file to be read
#' @param fileout Path and name of the csv file to be written
#' @param arr An array to be written to fileout
#' @param append Should an existing file be overwritten? By default, no.
#' @return A named list of arrays
#' 
#' @export
#' 
gv_readcsv <- function (filein)
{
	# first pass
	arr_nrows <- c()
	arr_names <- c()
	csv <- file(filein, 'r')
	line<- readLines(csv, n=1)
	while (length(line))
	{
		if (length(grep('^#', line, perl=T)))
		{
			words <- strsplit(line, ',\\s', perl=T)[[1]]
			arr_names <- append(arr_names, strsplit(words[1], '\\s', perl=T)[[1]][2])
			arr_nrows <- append(arr_nrows, as.numeric(words[2]))
		}
		line<- readLines(csv, n=1)
	}
	close(csv)
	# second pass
	arr <- list()
	n_skip <- 0
	n_rows <- cumsum(arr_nrows)
	for (i in 1:length(arr_nrows))
	{
		if (i != 1) { n_skip <- cumsum(n_rows[i-1]) }
		arr[[i]] <- as.matrix(
			read.csv(filein, header=F, skip=n_skip+i, nrows=arr_nrows[i])
		)
	}
	names(arr) <- arr_names
	return(arr)
}

#' Read and write a custom formatted CSV file
#' 
#' @inheritParams gv_sci_pix2vis_unsorted
#' 
#' @export
#' 
gv_writecsv <- function (fileout, arrname, arr, append=F)
{
  sz <- dim(arr)
  if (length(sz) == 0) { sz <- c(1,length(arr)); arr <- rbind(arr) }
  hdr <- sprintf('# %s, %d, %d', arrname, sz[1], sz[2])
  write.table(hdr, file=fileout, append=append,
              sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(arr, file=fileout, append=T,
              sep=',', row.names=FALSE, col.names=FALSE, quote=FALSE)
}
