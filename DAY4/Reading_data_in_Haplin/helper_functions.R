# Some updated functions - will be published on CRAN soon

f.get.which.gen.el <- function( cols, ncols.per.chunk ){
	if( !is.numeric( cols ) ){
		stop( "The 'cols' argument should be numeric!" )
	}

	which.gen.chunk <- vector( length( cols ), mode = "integer" )
	which.cols.chunk <- vector( length( cols ), mode = "integer" )
	for( i in 1:length( cols ) ){
		which.gen.chunk[ i ] <- ( cols[ i ] - 1 ) %/% ncols.per.chunk + 1
		which.cols.chunk[ i ] <- ( cols[ i ] - ( which.gen.chunk[ i ] - 1 )*ncols.per.chunk )%%( ncols.per.chunk + 1 )
	}
	
	return( list( chunk.no = which.gen.chunk, col.no = which.cols.chunk ) )
}
f.get.which.gen.el.names <- function( cols, gen.data ){
	if( !is.character( cols ) ){
		stop( "The 'cols' argument should be a character vector!" )
	}

	out.list <- sapply( cols, function( col.name ){
		for( chunk.no in 1:length( gen.data )){
			which.col <- match( col.name, colnames( gen.data[ chunk.no ] ), nomatch = -1 )
			if( which.col != -1 ){
				return( list( chunk.no = chunk.no, col.no = which.col ) )
			}
		}
	} )
	
	out.df <- data.frame( chunk.no = as.numeric( out.list[ "chunk.no", ] ), col.no = as.numeric( out.list[ "col.no", ] ) )
	return( out.df )
}

f.get.gen.data.cols <- function( gen.data, cols, by.colname = FALSE ){
	ncols.per.chunk <- ncol( gen.data[[ 1 ]] )
	if( by.colname ){
		gen.list.info <- f.get.which.gen.el.names( cols, gen.data )
	} else {
		gen.list.info <- f.get.which.gen.el( cols, ncols.per.chunk )
	}
	which.gen.chunk <- gen.list.info$chunk.no
	which.cols.chunk <- gen.list.info$col.no
	
	new.dim <- c( nrow( gen.data[[ 1 ]] ), length( cols ) )

	all.levels <- c()
	for( i in unique( which.gen.chunk ) ){
		all.levels <- union( all.levels, levels( gen.data[[ i ]] ) )
	}
	
	new.colnames <- c()
	gen.data.out <- ff::ff( NA, levels = all.levels, dim = new.dim, vmode = ff::vmode( gen.data[[ 1 ]] ) )
	for( i in 1:length( cols ) ){
		gen.data.out[ ,i ] <- gen.data[[ which.gen.chunk[ i ] ]][ ,which.cols.chunk[ i ] ]
		new.colnames <- c( new.colnames, 
			colnames( gen.data[[ which.gen.chunk[ i ] ]] )[ which.cols.chunk[ i ] ] )
	}
	colnames( gen.data.out ) <- new.colnames
	return( gen.data.out )
}

f.sel.markers <- function(n.vars, markers, family, split, ncols){
## SET NUMBER OF COLUMNS PER FAMILY (OR PER CHILD)
if(family == "mfc"){
	if(split) .t <- 6
	else .t <- 3
}else if(family == "c"){
	if(split) .t <- 2
	else .t <- 1
} else stop("Problem with family argument!", call. = F)
.test <- (ncols - n.vars) / .t 
.err <- F
#
if(!is.numeric(markers)){
	if(round(.test) != .test) .err <- T
	.markers <- 1:.test
	.sel <- 1:ncols
}else{
	## COMPUTE COLUMN NUMBERS OF MARKERS SPECIFIED IN markers ARGUMENT 
	.sel <- c(seq(length.out = n.vars), n.vars + as.numeric(t(outer((markers-1)*.t, 0:(.t - 1), "+")) + 1))
	if(max(.sel) > ncols) .err <- T
	.markers <- markers
}
if(.err){
	stop('It appears that the number of columns in the data file is wrong! \nCheck the datafile, the "n.vars", "sep" and "markers" arguments, etc....', call. = F)
}
#
##
attr(.sel, "markers") <- .markers
attr(.sel, "nloci") <- length(.markers)
return(.sel)
}

f.split.vector <- function(vec, split, na.strings = "NA"){
##
## EFFECTIVE SPLIT OF GENOTYPE VECTOR vec INTO TWO COLUMN MATRIX, SPLIT BY SEPARATOR split
## CHECKS THAT THERE'S ONLY ONE SPLIT. C;NA ETC IS SPLIT TO NA NA
##
#
## FIND ALL UNIQUE GENETIC VALUES, NA FIRST
.unique <- sort(unique(vec), na.last = F)
#
## TO SAVE TIME, SPLIT ONLY THE UNIQUE VALUES
.unique.split <- strsplit(.unique, split = split, fixed = T)
## JUST TO MAKE SURE, CHECK MISSING
.una <- which(is.na(.unique.split)) # SHOULD BE 1 IF MISSING, OTHERWISE EMPTY
if(length(.una) > 1) stop("Something's wrong with the data reading!\n")
#
## REPLACE THE SINGLE MISSING WITH TWO MISSING
if(length(.una) == 1) .unique.split[[.una]] <- c(NA_character_, NA_character_)
#
## CHECK THAT SPLITTING IS OK BY CHECKING THAT ALL GENETIC VALUES HAVE BEEN SPLIT INTO TWO
.ind <- sapply(.unique.split, length)
.ind <- which(.ind != 2)
if(length(.ind) > 0){
	cat("\n-----------------------\nSomething's wrong with the separator\nin the following data values:\n")
	print(.unique[.ind])
	stop('Correct data vector or check "split" argument!', call. = F)
}
#
## CONVERT TO MATRIX WITH TWO COLUMNS, ONE FOR EACH ALLELE
.unique.split <- t(sapply(.unique.split, function(x)x))
if(dim(.unique.split)[1] != length(.unique)) stop()
#
## MAKE SURE THERE ARE NO LEFTOVER MISSING OF THE TYPE NA;T ETC,
## RECODE THE T (OR WHATEVER) TO MISSING.
## THOSE OF TYPE NA;NA SHOULD ALREADY HAVE BEEN TAKEN CARE OF
.unique.split[.unique.split == na.strings] <- NA
.rna <- rowSums(is.na(.unique.split))
.unique.split[.rna > 0, ] <- NA
#
## EXPAND .unique.split TO MATCH THE vec LENGTH
.mat <- .unique.split[match(vec, .unique), ]
#
##
return(.mat)
}

f.split.matrix <- function(mat, split, tag.sep = "_"){
## DIMENSION AND DIMNAMES OF ORIGINAL MATRIX
.d <- dim(mat)
.dimnames <- dimnames(mat)
#
## RESHAPE TO CHARACTER VECTOR
.mat <- as.vector(mat)
#
## SPLIT ALLELES, RESULT IS A TWO-COLUMN MATRIX
.mat <- f.split.vector(.mat, split = split)
#
## SPLIT INTO TWO FILES, ONE FOR EACH ALLELE
.mat1 <- matrix(.mat[,1], nrow = .d[1], ncol = .d[2])
.mat2 <- matrix(.mat[,2], nrow = .d[1], ncol = .d[2])
#
## DIMENSION FOR JOINED DATA SET
.d <- dim(.mat1) * c(1,2)
## NEW JOINED DATA SET
.ut <- matrix(NA_character_, nrow = .d[1], ncol = .d[2])
.ut[, seq(1, .d[2], 2)] <- .mat1
.ut[, seq(2, .d[2], 2)] <- .mat2
#
## ADD CORRECT NAMES
row.names(.ut) <- rownames(mat)
.colnames <- outer(colnames(mat), 1:2, paste, sep = tag.sep)
.colnames <- as.vector(t(.colnames))
colnames(.ut) <- .colnames
#
##
return(.ut)
}

create.missingness.matrix <- function( 
	data.in = stop( "No data given!", call. = FALSE ),
	new.ids = stop( "No IDs given", call. = FALSE )
){

	id <- new.ids$ids
	pedIndex <- new.ids$pedIndex

	check.rows.only.NAs <- function( gen.data, ind ){
		sum.not.na.list <- lapply( gen.data, function( gen.el.ff ){
			gen.el <- gen.el.ff[ ind, ]
			which.na.level <- which( is.na( levels( gen.el ) ) )
			apply( gen.el, 1, function( row ){
				sum( as.numeric( row ) != which.na.level )
			} )
		})
		if( length( sum.not.na.list ) == 1 ){
			sum.not.na <- matrix( sum.not.na.list[[1]], ncol = 1 )
		} else {
			sum.not.na <- Reduce( cbind, sum.not.na.list )
			colnames( sum.not.na ) <- NULL
		}
		# this tells us which rows have only NAs
		rows.only.na <- rowSums( sum.not.na ) == 0
		return( rows.only.na )
	}

	# check the fathers first - which have genetic data?
	d.f <- match( pedIndex[ ,'id.father' ], id )
	d.f.NAs <- check.rows.only.NAs( data.in$gen.data, ind = d.f )
	ids.f.gen.data <- (id[ d.f ])[ !d.f.NAs ]
	
	# now, for the mothers...
	d.m <- match( pedIndex[ ,'id.mother' ], id )
	d.m.NAs <- check.rows.only.NAs( data.in$gen.data, ind = d.m )
	ids.m.gen.data <- (id[ d.m ])[ !d.m.NAs ]
	
	# and the children
	d.c <- match( pedIndex[ ,'id.child' ], id )
	d.c.NAs <- check.rows.only.NAs( data.in$gen.data, ind = d.c )
	ids.c.gen.data <- (id[ d.c ])[ !d.c.NAs ]
	
	# create a matrix where each member of the family has TRUE or FALSE
	#  (whether they have any genetic data)
	which.gen.data.fam <- as.data.frame( pedIndex, stringsAsFactors = FALSE )
	which.gen.data.fam$id.child <- FALSE
	which.gen.data.fam$id.father <- FALSE
	which.gen.data.fam$id.mother <- FALSE

	which.gen.data.fam$id.mother[ 
		match( ids.m.gen.data, pedIndex[ , "id.mother" ] ) ] <- TRUE
	which.gen.data.fam$id.father[ 
		match( ids.f.gen.data, pedIndex[ , "id.father" ] ) ] <- TRUE
	which.gen.data.fam$id.child[ 
		match( ids.c.gen.data, pedIndex[ , "id.child" ] ) ] <- TRUE

	return( which.gen.data.fam )
}

getFullTriads <- function( data.in = stop( "No data given!", call. = FALSE ),
	file.out = "my_data_onlyTriads",
	dir.out = ".", overwrite = NULL ){

	if( !is( data.in, "haplin.data" ) ||
		!all( names( data.in ) == Haplin:::.haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	format <- data.in$aux$info$filespecs$format
	
	if( format == "ped" ){
		# first - get all the IDs and check family structure
		new.ids <- f.check.unique.ids( data.in$cov.data )
		id <- new.ids$ids
		pedIndex <- new.ids$pedIndex

  	# check if there are any families that have missing members in pedIndex
  	any.NAs.families <- apply(pedIndex, 1, function(row){
  	  any(is.na(row))
  	})
  	pedIndex <- new.ids$pedIndex <- pedIndex[!any.NAs.families, ]
	
		which.gen.data.fam <- create.missingness.matrix( data.in, new.ids )

		# here are IDs of the members where the full family information is available
		full.triads <- apply( which.gen.data.fam[ ,-1 ], 1, all )
		pedIndex.triads <- pedIndex[ full.triads, ]
		
		# check how many families found
		if( nrow( pedIndex.triads ) == 0 ){
			stop( "No full triads found!", call. = FALSE )
		}
		
		full.triads.rows <- sort( c( match( pedIndex.triads[ ,'id.father' ], id ),
													match( pedIndex.triads[ ,'id.mother' ], id ),
													match( pedIndex.triads[ ,'id.child' ], id ) ) )
		return( genDataGetPart( data.in, design = "triad",
										rows = full.triads.rows,
										file.out = file.out, dir.out = dir.out, overwrite = overwrite )
					)
	} else if( format == "haplin" ){
		stop( "Not implemented yet", call. = FALSE )
	} else {
		stop( paste( "Unrecognized format:", format ), call. = FALSE )
	}
}

f.check.unique.ids <- function( data.cov ){
	all.id.c <- table( data.cov[ ,"id.c" ] )
	if( any( all.id.c > 1 ) ){
		cat( "   Creating unique IDs for individuals...\n" )
		orig.cov.colnames <- colnames( data.cov )
		data.cov <- t( apply( data.cov, 1, function( x ){
			if( x[ 4 ] == 0 | x[ 3 ] == 0 ){
				new.ids <- c( paste( x[ 1 ], x[ 2 ], sep = "_" ), x[ 3:4 ] )
			} else {
				new.ids <- paste( x[ 1 ], x[ 2:4 ], sep = "_" )
			}
			return( c( x[ 1 ], new.ids, x[ 5:length( x ) ] ) )
		} ) )
		colnames( data.cov ) <- orig.cov.colnames
		cat( "   ...done.\n" )
	}
	
	id <- data.cov[ ,"id.c" ]
	# sort the families and check coding
	pedIndex <- f.prep.pedIndex( data.cov )

	return( list( ids = id, pedIndex = pedIndex ) )
}

f.prep.pedIndex <- function( data.cov ){
## EXTRACT PED-INFORMATION

## EXTRACT FAMILY INFORMATION FROM PED FILE
# cat("Extracting family, sex and case/control information from ped file...\n")
.fam <- as.data.frame( data.cov, stringsAsFactors = FALSE )

names( .fam ) <- c( "family", "id", "father", "mother", "sex", "cc" )
.pheno <- as.matrix( .fam[ , c( "id", "sex", "cc" ) ], mode = "character" )
.fam <- .fam[, c( "family", "id", "father", "mother" ) ]

.sex.u <- sort( unique( .pheno[ ,"sex" ] ) )
if( length( .sex.u ) > 2 ){
	stop( "More than 2 different codes in the sex column (column 5)", call. = F )
}
if( any( !is.element( .sex.u, c("0", "1") ) ) ){
	if( all( is.element( .sex.u, c("1", "2") ) ) ){
		.pheno[ .pheno[ ,"sex" ] == "2","sex" ] <- "0" # RECODE FEMALES
	}else{
		stop("Invalid codes in the sex column (column 5)", call. = F)
	}
}
.pheno <- as.data.frame( .pheno, stringsAsFactors = FALSE )

## CREATE INDEXING TO BE USED LATER WHEN CONVERTING TO HAPLIN FORMAT
.pedIndex <- f.make.index( .fam, output = "ids" )

# return( invisible() )
return( .pedIndex )
}

f.make.index <- function(vardata, output = "line numbers"){
## CHECK family VARIABLE
cat("\nChecking family and id variables...\n")
if( any( table( vardata$family, exclude = NULL ) > 3) ){
	warning("Found family size larger than 3!  Will extract trios from general pedigree.", call. = F)
}

if( .test <- ( any( vardata$family == "0") | any( vardata$id == "0" ) ) ){
	stop( paste( 'Cannot use "0" in family or id code!\n', 'Found on line(s): ', paste( which( .test ), collapse = " " ), sep = "" ), call. = F )
}

if( .test <- ( any( is.na( vardata$family ) ) | any( is.na( vardata$id ) ) ) ){
	stop( paste( 'Cannot have missing values in family or id variable\n', 'Found on line(s): ', paste( which( .test ), collapse = " "), sep = "" ), call. = F )
}

## RECODE ZEROES TO MISSING
vardata$mother[vardata$mother == "0"] <- NA
vardata$father[vardata$father == "0"] <- NA
#
## CREATE A VARIABLE THAT UNIQUELY IDENTIFIES INDIVIDUALS
.tagit <- "<>" # I BET THAT ONE'S UNIQUE!
.tag <- paste(vardata$family, vardata$id, sep = .tagit)
.tag.mother <- paste(vardata$family, vardata$mother, sep = .tagit)
.tag.father <- paste(vardata$family, vardata$father, sep = .tagit)
## RETAIN MISSING
.tag.mother[is.na(vardata$mother)] <- NA
.tag.father[is.na(vardata$father)] <- NA
#
## CHECK FOR DUPLICATES
.dupl <- duplicated(.tag)
if(any(.dupl)) {
	cat("\n")
	.mess <- paste("Individual id appears several times within one family!\nFor instance, family ", vardata$family[.dupl][1], " contains more than one of individual ", vardata$id[.dupl][1], ".", sep = "")
	stop(.mess, call. = F)
}
#
## IDENTIFY MOTHERS AND FATHERS
.is.mother <- is.element(.tag, .tag.mother)
.is.father <- is.element(.tag, .tag.father)
#
## ELIMINATE PARENT CODINGS THAT REFER TO NON-EXISTENT INDIVIDUALS,
## SAY, IF LINES OF THE FILE HAVE BEEN REMOVED
.rem.mother <- !is.element(.tag.mother, c(NA, .tag))
.rem.father <- !is.element(.tag.father, c(NA, .tag))
.sm <- sum(.rem.mother)
.sf <- sum(.rem.father)
if(.sm + .sf > 0){
	.tag.mother[.rem.mother] <- .tag.father[.rem.father] <- NA
	.mess <- paste(.sm, " mother code(s) and ", .sf, " father code(s) refer to non-existing individuals and have been set to missing", sep = "")
	warning(.mess, call. = F)
}
#
## DEFINE AND IDENTIFY "CHILD"
## 1) Child = person having a mother or father
.is.child.1 <- !is.na(.tag.mother) | !is.na(.tag.father)
## 2) or = person not having mother nor father but also does not have any children, i.e. "all alone"
.is.child.2 <- !.is.child.1 & (!.is.mother & !.is.father)
#
.is.child <- .is.child.1 | .is.child.2
#  
## SOME CHECKING
if(any(.test <- .is.father & .is.mother)){
	stop(paste("Sorry, the same individual cannot be both father and mother!\n", "Found on line(s): ", paste(which(.test), collapse = " "), sep = ""), call. = F)
}
#
## FIND LINE NUMBER OF CHILD AND OF ITS PARENTS
.line.child <- which(.is.child)
.line.mother <- match(.tag.mother, .tag)[.is.child]
.line.father <- match(.tag.father, .tag)[.is.child]
#
##
if(output == "line numbers"){
	## OUTPUT GIVES LINE NUMBERS, NOT IDs.
	## THIS COULD BE USEFUL FOR GenABEL SINCE FAMILY IDS ARE NOT RETAINED
	## (BUT SHOULD MAKE SURE SUBSETTING IN GenABEL IS TAKEN INTO ACCOUNT WHEN CONVERTING TO HAPLIN)
	#
	## JOIN LINE NUMBERS INTO MATRIX
	.ut <- cbind(line.child = .line.child, line.mother = .line.mother, line.father = .line.father)
}else if(output == "ids"){
	## USES SAME IDS AS IN FILE. REQUIRES IDS TO BE UNIQUE!
	## CAN BE USED WITH GenABEL SINCE THEY INSIST THAT IDS SHOULD BE UNIQUE.
	if( any( duplicated( vardata$id ) ) ){
		stop("Duplicated individual id", call. = F)
	}
	.fam.child <- vardata$family[.line.child]
	.fam.mother <- vardata$family[.line.mother]
	.fam.father <- vardata$family[.line.father]
	#
	.id.child <- vardata$id[.line.child]
	.id.mother <- vardata$id[.line.mother]
	.id.father <- vardata$id[.line.father]
	#
	if(any(.fam.child != .fam.mother, na.rm = T) | any(.fam.child != .fam.father, na.rm = T)){
		stop( "Problem with family identification!", call. = F )
	}

	## JOIN FAMILY ID, CHILD ID, MOTHER ID, FATHER ID INTO MATRIX
	.ut <- cbind(family = .fam.child, id.child = .id.child, id.mother = .id.mother, id.father = .id.father)
}else if(output == "tags"){
	## USES TAGS CREATED FROM FAMILY ID COMBINED WITH IDS AS IN FILE.
	## CAN BE USED IN GENERAL HAPLIN CONVERSIONS WHICH DO NOT REQUIRE ID ITSELF TO BE UNIQUE
	## BUT WILL NOT WORK WITH GenABEL SINCE LATTER DOES NOT RETAIN FAM IDS
	.fam.child <- vardata$family[.line.child]
	.fam.mother <- vardata$family[.line.mother]
	.fam.father <- vardata$family[.line.father]
	#
	.id.child <- .tag[.line.child]
	.id.mother <- .tag[.line.mother]
	.id.father <- .tag[.line.father]
	#
	if(any(.fam.child != .fam.mother, na.rm = T) | any(.fam.child != .fam.father, na.rm = T)) stop("Problem with family identification!", call. = F)
	#
	## JOIN FAMILY ID, CHILD ID, MOTHER ID, FATHER ID INTO MATRIX
	.ut <- cbind(family = .fam.child, id.child = .id.child, id.mother = .id.mother, id.father = .id.father)
}
#
##
return(.ut)
}

f.make.out.filename <- function( file.in, file.out, dir.out, root = "gen", overwrite = NULL ){
	if( is.null( file.out ) ){
		file.in.base <- basename( file.in )
		fname.vec <- unlist( strsplit( file.in.base, split = ".", fixed = TRUE ) )
		l.fname.vec <- length( fname.vec )
		if( l.fname.vec > 1 ){
			file.out <- paste( fname.vec[ -l.fname.vec ], collapse = "_" )
		} else {
			file.out <- file.in.base
		}
	}
	
	file.out.base <- paste0( file.out, "_", root )
	file.out.ff <- file.path( dir.out, paste0( file.out.base, ".ffData") )
	file.out.aux <- file.path( dir.out, paste0( file.out.base, ".RData" ) )
	if( file.exists( file.out.ff ) | file.exists( file.out.aux ) ){
		cat( "The output file(s) exist! \n" )
		if( is.null( overwrite ) ){
			answer <- readline( paste( "Do you want to overwrite file(s)? (y/n)" ) )
			if( tolower( answer ) != "y" ){
				stop( "Stopped without overwriting files.", call. = FALSE )
				return( NULL )
			}
		} else if( !is.logical( overwrite ) | !overwrite ){
			stop( "Stopped without overwriting files.", call. = FALSE )
			return( NULL )
		} else if( overwrite ){
			ff::ffdrop( file.path( dir.out, file.out.base ) )
		}
	}
	
	return( list( file.out.base = file.out.base, file.out.ff = file.out.ff, file.out.aux = file.out.aux ) )
}

f.create.snp.names <- function( map.file, ncol, format, design ){
	cat( "Reading the marker names... \n" )

	if( format == "ped" ){
		ncol.per.locus <- 2
	} else if( design %in% c( "triad", "cc.triad" ) ){
		ncol.per.locus <- 6
	} else {
		ncol.per.locus <- 2
	}

	marker.names <- c()
	if( !is.null( map.file ) ){
		marker.names <- read.table( map.file, header = TRUE, stringsAsFactors = FALSE )
		if( ( nrow( marker.names ) * ncol.per.locus ) != ncol ){
			marker.names <- c()
		}
	}
	if( length( marker.names ) == 0 ){
		warning( "No map file given, map file empty or the number of map file rows not equal to the number of markers in data; will generate dummy marker names.", call. = FALSE )
		marker.names <- paste( "m", 1:( ncol/ncol.per.locus ), sep = "" )
	} else {
		marker.names <- marker.names[ ,2 ]
	}

	# all the marker names will start with l_
	gen.data.colnames <- paste( "l", as.character( marker.names ), sep = "_" )
	markers1 <- paste( gen.data.colnames, "a", sep = "_" )
	markers2 <- paste( gen.data.colnames, "b", sep = "_" )

	if( format == "ped" | ( format == "haplin" & design == "cc" ) ){
		gen.data.colnames <- as.vector( rbind( markers1, markers2 ) )
	} else if( format == "haplin" & design %in% c( "triad", "cc.triad" ) ){
		labs <- c("m", "f", "c")
		marker.names.a <- as.vector( t( outer( markers1, labs, paste, sep = "_" ) ) )
		marker.names.b <- as.vector( t( outer( markers2, labs, paste, sep = "_" ) ) )
		gen.data.colnames <- as.vector( rbind( marker.names.a, marker.names.b ) )
	} else {
		stop( "Problem with design and format!", call. = FALSE )
	}
	cat( "...done\n" )
	
	return( list( gen.data.colnames = gen.data.colnames, marker.names = marker.names ) )
}

genDataGetPart <- function( data.in = stop( "No data given!", call. = FALSE ), design = stop( "Design type must be given!" ), markers, indiv.ids, rows, cc, sex, file.out = "my_data_part", dir.out = ".", overwrite = NULL, ... ){
	#---- checking the input params ---------------------
	files.list <- f.make.out.filename( file.out = file.out, dir.out = dir.out, overwrite = overwrite )
	
	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == Haplin:::.haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	design.list <- get( ".design.list", envir = Haplin:::.haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}
	#---- done checking the input params ----------------
	#----------------------------------------------------
	# get all arguments
	all.args <- match.call()[-1]
	all.args.names <- names( all.args )
	# arguments that do not define selection of data subset:
	non.sel.args <- c( "data.in", "file.out", "dir.out", "overwrite", "design" )
	
	which.non.sel.args.given <- c()
	for( i in 1:length( non.sel.args ) ){
		which.cur.arg <- match( non.sel.args[ i ], all.args.names )
		if( !is.na( which.cur.arg ) ){
			which.non.sel.args.given <- c( which.non.sel.args.given, which.cur.arg )
		}
	}
	# arguments that define selection of data subset:
	selection.args <- all.args[ -which.non.sel.args.given ]

	# this is the maximum set of arguments for subset selection, based on the read in covariate data
	format <- data.in$aux$info$filespecs$format
# 	if( format == "ped" ){
# 		correct.sel.args <- union( names( formals() ), colnames( data.in$cov.data )[ -(1:4) ] )
# 	} else { # haplin file
		correct.sel.args <- union( names( formals() ), colnames( data.in$cov.data ) )
# 	}
	correct.sel.args <- correct.sel.args[ -( match( c( non.sel.args, "..." ),
		correct.sel.args ) ) ]
	
	max.rows <- nrow( data.in$gen.data[[1]] )
	all.rows <- 1:max.rows
	all.markers <- sum( unlist( lapply( data.in$gen.data, ncol ) ) )
	all.cols <- 1:all.markers
	
	subset.rows <- all.rows
	subset.cols <- all.cols
	is.subset.cols <- FALSE
	is.subset.rows <- FALSE
	family <- "c"
	if( design %in% c( "triad", "cc.triad" ) & format == "haplin" ){
		family <- "mfc"
	}
	
	cat( "Provided arguments:\n" )
	for( arg in names( selection.args ) ){
		if( arg == "markers" ){
			val <- eval( all.args$markers, envir = rlang::caller_env() )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- chosen markers:", print.val, "\n" )
			if( any( val > all.markers ) | any( val < 0 ) ){
				stop( "wrong markers chosen; not in range of given data: max ", all.markers, "!\n", call. = FALSE )
			}
			if( identical( val, all.cols ) ){
				warning( "this selection is equal to choosing all markers!", call. = FALSE )
			} else {
				subset.cols <- f.sel.markers( n.vars = 0, markers = val, family = family, split = TRUE, ncols = all.markers )
				is.subset.cols <- TRUE
			}
			
		} else if( arg == "indiv.ids" ){
			val <- eval( all.args$indiv.ids, envir = rlang::caller_env() )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- individual IDs:", print.val, "\n" )
			
			if( is.null( data.in$cov.data ) |
				!( "id.c" %in% colnames( data.in$cov.data ) ) ){
				# this can happen when the haplin-formatted file had no extra info
				stop( "There is no information on individual IDs in the data!", call. = FALSE )
			}
			
			chosen.indiv.logic <- val %in% as.character( data.in$cov.data[ ,"id.c" ] )
			if( !any( chosen.indiv.logic ) ){
				stop( "wrong individual IDs chosen!", call. = FALSE )
			} else if( sum( chosen.indiv.logic ) == max.rows ){
				warning( "this selection is equal to choosing all the individuals available in dataset!", call. = FALSE )
			} else {
				chosen.indiv <- do.call( c, sapply( val, function( x ){
					list( which( as.character( data.in$cov.data[ ,"id.c" ] ) == x ) )
				} ) )
				subset.rows <- intersect( subset.rows, chosen.indiv[ !( is.na( chosen.indiv ) ) ] )
				is.subset.rows <- TRUE
			}
			
		} else if( arg == "rows" ){
			val <- eval( all.args$rows, envir = rlang::caller_env() )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- rows chosen:", print.val, "\n" )
			if( any( val < 1 ) | any( val > nrow( data.in$cov.data ) ) ){
				stop( "wrong rows, not in range of given data: max ", max.rows, "\n", call. = FALSE )
			}
			if( identical( val, all.rows ) ){
				warning( "this selection is equal to choosing all the rows available in dataset!", call. = FALSE )
			} else {
				subset.rows <- intersect( subset.rows, val )
				is.subset.rows <- TRUE
			}
			
		} else if( !( arg %in% correct.sel.args ) ){
 			stop( "the given argument: ", arg, " is not recognizable!", call. = FALSE )
 			
		} else { #any other argument if additional covariate data are loaded
			val <- eval( all.args[[ match( arg, names( all.args ) ) ]], envir = rlang::caller_env() )
			print.val <- paste( val, collapse = ", " )
			if( length( val ) > 20 ){
				print.val <- paste( paste( head( val ), collapse = ", "), "...", paste( tail( val ), collapse = ", " ) )
			}
			cat( " --- argument: ", arg,", chosen values: ", print.val, "\n", sep = "" )
			avail.vals <- unique( as.numeric( data.in$cov.data[ ,arg ] ) )
			if( !all( val %in% avail.vals ) ){
				stop( "wrong argument(", arg, ") chosen, available: ", paste( avail.vals, sep = " ", collapse = " " ), "!", call. = FALSE )
			}
			if( all( avail.vals %in% val ) ){
				warning( "this selection is equal to choosing all the possible values of the given covariate (", arg, ")!", call. = FALSE )
			} else {
				chosen.indiv <- do.call( c, sapply( val, function( x ){
					list( which( as.numeric( data.in$cov.data[ ,arg ] ) == x ) )
				} ) )
				subset.rows <- intersect( subset.rows, chosen.indiv )
				is.subset.rows <- TRUE
			}
			
	 	}
	}
	if( is.subset.rows | is.subset.cols ){
		cat( "INFO: Will select ", length( subset.rows ), " rows and ", length( subset.cols ), " columns.\n", sep = "" )
	} else {
		stop( "No subset selected.\n", call. = FALSE )
	}
	
	if( is.subset.cols ){
		# check how many chunks will be needed
		nb.cols.per.chunk <- get( ".nb.cols.per.chunk", envir = Haplin:::.haplinEnv )
		nb.chunks <- ceiling( length( subset.cols ) / nb.cols.per.chunk )
		
		gen.data.col.wise <- sapply( 1:nb.chunks, function( x ){
			first.col <- ( nb.cols.per.chunk*( x - 1 ) + 1 )
			last.col <- min( ( nb.cols.per.chunk*x ), length( subset.cols ) )
			cur.cols <- subset.cols[ first.col:last.col ]
			list( f.get.gen.data.cols( data.in$gen.data, cur.cols ) )
		} )
	} else {
		gen.data.col.wise <- data.in$gen.data
	}
	
	cov.data.in <- NULL
	if( is.subset.rows ){
		# need to choose from both gen.data and cov.data
		gen.data.col.wise <- lapply( gen.data.col.wise, function( x ){
			sub <- x[ subset.rows, ]
			out <- ff::ff( sub, levels = levels( sub ), dim = dim( sub ), vmode = ff::vmode( x ) )
			colnames( out ) <- colnames( sub )
			return( out )
		})
		
		if( !is.null( data.in$cov.data ) ){
			cov.data.in <- data.in$cov.data[ subset.rows, ]
		}
	} else if( !is.null( data.in$cov.data ) ){
		cov.data.in <- data.in$cov.data
	}
	data.out <- list( cov.data = cov.data.in, gen.data = gen.data.col.wise, aux = data.in$aux )
	class( data.out ) <- class( data.in )
	
	# update the marker names
	if("markers" %in% names(selection.args)){
	  which.markers.retain <- eval(all.args$markers, envir = rlang::caller_env())
	  data.out$aux$marker.names <- data.in$aux$marker.names[which.markers.retain]
	}
	
	## saving the chosen part of the data
	cat( "Saving data... \n" )
	cur.names <- c()
	for( i in 1:length( gen.data.col.wise ) ){
		cur.name <- paste( get( ".gen.cols.name", envir = Haplin:::.haplinEnv ), i, sep = "." )
		assign( cur.name, gen.data.col.wise[[i]] )
		cur.names <- c( cur.names, cur.name )
	}
	aux <- data.out$aux
	save.list <- c( cur.names, "aux" )
	if( !is.null( cov.data.in ) ){
		save.list <- c( save.list, "cov.data.in" )
	}
	ff::ffsave( list = save.list, file = file.path( dir.out, files.list$file.out.base ) )
	cat( "... saved to files: ", files.list$file.out.ff, ", ", files.list$file.out.aux, "\n", sep = "" )

	return( data.out )
}

genDataRead <- function( file.in = stop( "Filename must be given!", call. = FALSE ), file.out = NULL, dir.out = ".", format = stop( "Format parameter is required!" ), header = FALSE, n.vars, cov.file.in, cov.header, map.file, allele.sep = ";", na.strings = "NA", col.sep = "", overwrite = NULL ){
	## checking the input arguments
	if( !file.exists( file.in ) ){
		stop( "The given file (", file.in, ") doesn't exist! Check and try again.", call. = FALSE )
	}
	if( !missing( cov.file.in ) ){
		if( !file.exists( cov.file.in ) ){
			stop( "The given file (", cov.file.in, ") doesn't exist! Check and try again.", call. = FALSE )
		}
	}
	if( !missing( map.file ) ){
		if( !file.exists( map.file ) ){
			stop( "The given map.file (", map.file, ") doesn't exist! Check and try again.", call. = FALSE )
		}
	} else {
	  map.file <- NULL
	}

	files.list <- f.make.out.filename( file.in, file.out, dir.out = dir.out, overwrite = overwrite )
	
	if( !( format %in% c("haplin", "ped") ) ){
		stop( "Given format (", format,") not recognized! Format has to be one of: ", paste( c("haplin", "ped"), collapse = ", " ), call. = FALSE )
	}
	
	## some special checks if the input file is in haplin format:
	if( format == "haplin" ){
		if( ( na.strings == col.sep ) | ( ( na.strings != "" ) & ( na.strings == allele.sep ) ) ){
			stop( 'The "na.strings" argument should be different from "col.sep" and "allele.sep" arguments!', call. = F )
		}
		
		split <- F
		if( ( col.sep != "" & col.sep == allele.sep ) | allele.sep == " " ){
			split <- T
		}
	}
	
	if( !missing( cov.header ) ){
		if( !is.character( cov.header ) ){
			stop( "'cov.header' is specified, but it's not a character vector!", call. = FALSE )
		}
	}
	
	skip.first <- 0
	if( header ){
		header.line <- scan( file = file.in, what = "character", nlines = 1, strip.white = TRUE, sep = col.sep, na.strings = na.strings )
		skip.first <- 1
	}
	
	## read the first line of the file and check the number of columns
	first.line <- scan( file = file.in, what = "character", nlines = 1, strip.white = TRUE, sep = col.sep, na.strings = na.strings, skip = skip.first )
	nb.cols <- length( first.line )

	if( missing( n.vars ) ){
		if( format == "haplin" ){
			stop( 'The value of "n.vars", i.e. the number\n of covariate columns ahead of genetic data,\n should be specified when format = "haplin".\n Set n.vars = 0 if there are no covariates in haplin file.', call. = FALSE )
		}
		if( format == "ped" ){
			n.vars <- 6
		}
		cat( "'n.vars' was not given explicitly and will be set to ", n.vars, " based on the format given.\n", sep = "" )
	} else if( n.vars < 0 | n.vars >= nb.cols ){
		stop( "The 'n.vars' value (", n.vars, ") is not valid!", call. = FALSE )
	}

	if( format == "haplin" & n.vars > 0 & missing( cov.header ) ){
		cat( "The format of the file is 'haplin' with covariate data but no names of the covariate data is given. Will generate dummy names.\n" )
	}

	## the chunk size that would still fit in the memory
# 	nb.lines.per.chunk <- get( ".nb.lines.per.chunk", envir = Haplin:::.haplinEnv )
	nb.lines.per.chunk <- ceiling( 100000000 / nb.cols )

	## open the file
	in.file <- file( file.in, "r" )

	gen.data.in.ffdf <- list()
	cov.data.in <- c()
	gen.levels <- c()

	nb.rows.tot <- 0

	i <- 1
	cat( "Reading the data in chunks...\n" )
	cat( " -- chunk ", i, " --\n", sep = "" )
	cur.chunk <- matrix( scan( in.file, what = "character", nlines = nb.lines.per.chunk, sep = col.sep, na.strings = na.strings ), ncol = nb.cols, byrow = TRUE )

	## reading in chunks and creating ff object for each chunk
	while( length( cur.chunk ) != 0 ){
		if( n.vars > 0 ){
			cov.data.in <- rbind( cov.data.in, cur.chunk[ ,1:n.vars, drop = F ] )
			cur.chunk <- cur.chunk[ ,-(1:n.vars), drop = F ]
		}

		if( format == "haplin"){
			## SPLIT ALLELES IF NOT ALREADY DONE. IF AT LEAST ONE ALLELE IS MISSING, SET THE OTHER ONE TO MISSING, TOO
			if( !split ){
				cur.chunk <- f.split.matrix( cur.chunk, split = allele.sep )
			}else{
			## IF SPLIT ALREADY (MAY NOT BE NECESSARY TO SPLIT UP, BUT THIS CAME BEFORE f.split.matrix...)
				## KEEP DIMENSION
				d <- dim( cur.chunk )
				.data.gen1 <- cur.chunk[, seq(1, d[2], 2), drop = F]
				.data.gen2 <- cur.chunk[, seq(2, d[2], 2), drop = F]
				
				## SET TO MISSING ALL THOSE WITH AT LEAST ONE MISSING ALLELE
				.is.na <- is.na(.data.gen1) | is.na(.data.gen2)
				.data.gen1[.is.na] <- NA
				.data.gen2[.is.na] <- NA
				
				## PUT BACK TOGETHER AGAIN
				## DIMENSION FOR JOINED DATA SET
				.d <- dim(.data.gen1) * c(1,2)
				## NEW JOINED DATA SET
				cur.chunk <- matrix(NA_character_, nrow = d[1], ncol = d[2])
				cur.chunk[, seq(1, d[2], 2)] <- .data.gen1
				cur.chunk[, seq(2, d[2], 2)] <- .data.gen2
# 				row.names(.data.gen) <- 1:(dim(.data.gen))[1]
				
				rm(.data.gen1)
				rm(.data.gen2)
# 				gc()
			}
		}
	
		cur.levels <- unique( as.vector( cur.chunk ) )
		recode.na <- FALSE
		if( 0 %in% cur.levels ){
			recode.na <- TRUE
			na.symbol <- 0
		} else if( "0" %in% cur.levels ){
			recode.na <- TRUE
			na.symbol <- "0"
		}
		if( recode.na ){
			cur.chunk[ cur.chunk == na.symbol ] <- NA
			cur.levels[ cur.levels == na.symbol ] <- NA
		}
		gen.levels <- as.character( union( gen.levels, cur.levels ) )

		tmp.ff <- ff::as.ff( cur.chunk, vmode = Haplin:::.haplinEnv$.vmode.gen.data, levels = gen.levels )

		gen.data.in.ffdf <- c( gen.data.in.ffdf, list( tmp.ff ) )
		nb.rows.tot <- nb.rows.tot + nrow( tmp.ff )

		rm( cur.chunk, tmp.ff )
# 		gc()

		i <- i + 1
		cat( " -- chunk ", i, " -- \n", sep = "" )
		cur.chunk <- matrix( scan( in.file, what = "character", nlines = nb.lines.per.chunk ), ncol = nb.cols, byrow = TRUE )
	}
	cat( "... done reading.\n" )

	close( in.file )

	## re-organize - it's much better to have a list with different column-chunks
	nb.cols.per.chunk <- get( ".nb.cols.per.chunk", envir = Haplin:::.haplinEnv )
	nb.cols.gen.data <- ncol( gen.data.in.ffdf[[1]] )
	nb.col.chunks <- ceiling( nb.cols.gen.data / nb.cols.per.chunk )
	gen.list.length <- length( gen.data.in.ffdf )
	gen.data.col.wise <- list()
	
	design <- "cc"
	if( format == "haplin" ){
		design <- "triad"
	}
	tot.gen.ncol <- ncol( gen.data.in.ffdf[[ 1 ]] )
	gen.data.colnames <- f.create.snp.names( map.file, ncol = tot.gen.ncol, format = format, design = design )
	marker.names <- gen.data.colnames$marker.names
	gen.data.colnames <- gen.data.colnames$gen.data.colnames

	cat( "Preparing data...\n" )
	for( i in 1:nb.col.chunks ){
		cur.cols <- ( ( i-1 )*nb.cols.per.chunk + 1 ):( min( i*nb.cols.per.chunk, nb.cols.gen.data ) )
		tmp.gen.data <- ff::ff( vmode = Haplin:::.haplinEnv$.vmode.gen.data, levels = gen.levels, dim = c( nb.rows.tot, min( nb.cols.per.chunk, max( cur.cols ) - min( cur.cols ) + 1 ) ) )
		
		prev.rows <- 0
		for( j in 1:gen.list.length ){
			cur.rows <- ( prev.rows + 1 ):( prev.rows + nrow( gen.data.in.ffdf[[j]] ) )
			tmp.gen.data[ cur.rows, ] <- gen.data.in.ffdf[[j]][ ,cur.cols ]
			prev.rows <- max( cur.rows )
		}
		colnames( tmp.gen.data ) <- gen.data.colnames[ cur.cols ]
		
		gen.data.col.wise <- c( gen.data.col.wise, list( tmp.gen.data ) )
		rm( tmp.gen.data )
	}
	cat( "... done preparing\n" )

	rm( gen.data.in.ffdf )

	cov.data.colnames <- c()
	if( n.vars > 0 ){
		if( format == "ped" ){
			## lookup in the package environment
			cov.data.colnames <- get( ".cov.data.colnames", envir = Haplin:::.haplinEnv )
		} else if( !missing( header ) ) {
			cov.data.colnames <- header.line[ 1:n.vars ]
		} else {
			cov.data.colnames <- paste0( "cov.", 1:n.vars )
		}
	}
	
	## reading additional data (if given)
	cov.n.vars <- 0
	if( !missing( cov.file.in ) ){
		cat( "Reading covariate file... \n" )
		if( missing( cov.header ) ){
			cat( "    'cov.header' not given - assuming the first line is the header...\n" )
			cov.file.header <- TRUE
		} else {
			cov.file.header <- FALSE
		}

		cov.add.data <- read.table( cov.file.in, header = cov.file.header, stringsAsFactors = FALSE )
		if( missing( cov.header ) ){
			cov.header <- colnames( cov.add.data )
		}
		cov.data.colnames <- c( cov.data.colnames, cov.header )
		if( nrow( cov.add.data ) != nb.rows.tot ){
			stop( "The number of rows in the additional covariate data (", nrow( cov.add.data ), ") doesn't match the number of rows in the main file with genetic data (", nb.rows.tot, ")!", call. = FALSE )
		}
		
		if( !cov.file.header & ( ncol( cov.add.data ) != length( cov.header ) ) ){
			stop( "The length of given 'cov.header' names (", length( cov.header ), ") doesn't match the number of all the covariate data columns (", ncol( cov.add.data ), ")! Check and try again." )
		}
		
		if( length( cov.data.in ) != 0 ){
			cov.data.in <- cbind( cov.data.in, cov.add.data )
		} else {
			cov.data.in <- cov.add.data
		}
		cov.n.vars <- ncol( cov.add.data )
		cat( "...done\n" )
	} else if( !missing( cov.header ) ) {
		# if no additional file with covariates is read but there are some covariates in the
		# main file with genotypes and one wants to override the default values
		if( format == "haplin" ){
			cov.n.vars <- 0
			if( n.vars != length( cov.header ) ){
				stop( "The length of given 'cov.header' names (", length( cov.header ), ") doesn't match the number of all the covariate data columns (", n.vars, ")! Check and try again." )
			}else{
				cov.data.colnames <- cov.header
			}
		} else if( format == "ped" ){
			cov.n.vars <- n.vars - 6
			if( cov.n.vars != length( cov.header ) ){
				stop( "The length of given 'cov.header' names (", length( cov.header ), ") doesn't match the number of all the covariate data columns (", cov.n.vars, ")! Check and try again." )
			}else{
				cov.data.colnames[ 7:n.vars ] <- cov.header
			}
		}
	}
	n.vars <- n.vars + cov.n.vars
	
	## saving the data in the .RData and .ffData files
	cat( "Saving data...\n" )
	cur.names <- c()
	for( i in 1:length( gen.data.col.wise ) ){
		cur.name <- paste( get( ".gen.cols.name", envir = Haplin:::.haplinEnv ), i, sep = "." )
		assign( cur.name, gen.data.col.wise[[i]] )
		cur.names <- c( cur.names, cur.name )
	}
	save.list <- cur.names
	if( n.vars > 0 ){
		if( n.vars == 1 ){
			cov.data.in <- matrix( cov.data.in, ncol = 1 )
		}
		colnames( cov.data.in ) <- cov.data.colnames
		save.list <- c( save.list, "cov.data.in" )
	} else {
		cov.data.in <- NULL
	}

	## Include an info object
	.info <- list()
	.info$filename <- c( file.in = file.in )
	if( !missing( cov.file.in ) ){
		.info$filename <- c( .info$filename, cov.file.in = cov.file.in )
	}
	.info$filespecs$n.vars <- n.vars
	.info$filespecs$sep <- col.sep
	.info$filespecs$allele.sep <- allele.sep
	.info$filespecs$na.strings <- na.strings
	.info$filespecs$format <- format
	aux <- list( info = .info, class = "haplin.data" )
	aux$marker.names <- marker.names
	aux$map.filename <- map.file

	save.list <- c( save.list, "aux" )
	
	ff::ffsave( list = save.list, file = file.path( dir.out, files.list$file.out.base ) )
	cat( "... saved to files: ", files.list$file.out.ff, ", ", files.list$file.out.aux, "\n", sep = "" )

	data.out <- list( cov.data = cov.data.in, gen.data = gen.data.col.wise, aux = aux )
	class( data.out ) <- aux$class
	
	return( data.out )
}

genDataPreprocess <- function( data.in = stop( "You have to give the object to preprocess!" ), map.file, design = "triad", file.out = "data_preprocessed", dir.out = ".", ncpu = 1, overwrite = NULL ){
	# checking input parameters:
	if( !missing( map.file ) ){
		if( !file.exists( map.file ) ){
			stop( "It seems like the map.file (", map.file, ") doesn't exist! Check and try once more.", call. = FALSE )
		}
	} else if( !is.null( data.in$aux$map.filename ) ){
		map.file <- data.in$aux$map.filename
	} else {
		map.file <- NULL
	}
	
	.info <- data.in$aux$info
	.format <- .info$filespecs$format
	if( !( .format %in% c( "haplin", "ped" ) ) ){
		stop( "Given format (", format,") not recognized! Format has to be one of: ", paste( c("haplin", "ped"), collapse = ", " ), call. = FALSE )
	}
	
	design.list <- get( ".design.list", envir = Haplin:::.haplinEnv )
	if( !( design %in% design.list ) ){
		stop( "Given design(", design,") not recognized! Design has to be one of: ", paste( design.list, collapse = ", " ), call. = FALSE )
	}
	.info$model$design <- design
	
	files.list <- f.make.out.filename( file.out = file.out, dir.out = dir.out, overwrite = overwrite )
	

	if( !is( data.in, "haplin.data" ) ||
	  !all( names( data.in ) == Haplin:::.haplinEnv$.haplin.data.names ) ){
		stop( "The input data is not in the correct format!", call. = FALSE )
	}
	
	if( ncpu < 1 ){
		cat( " You set 'ncpu' to a number less than 1 - resetting it to 1.\n" )
		ncpu <- 1
	}
	#--------------------------------------------

	# creating SNP names (dummy names)
	tot.gen.ncol <- sum( unlist( lapply( data.in$gen.data, ncol ) ) )
	gen.data.colnames <- f.create.snp.names( map.file, tot.gen.ncol, format = .format, design = design )
	marker.names <- gen.data.colnames$marker.names
	gen.data.colnames <- gen.data.colnames$gen.data.colnames

	cur.col <- 1
	for( i in 1:length( data.in$gen.data ) ){
		cur.ncol <- ncol( data.in$gen.data[[ i ]] )
		colnames( data.in$gen.data[[ i ]] ) <- gen.data.colnames[ cur.col:( cur.col + cur.ncol - 1 ) ]
		cur.col <- cur.col + cur.ncol
	}

	# re-organizing data from PED format to internal haplin
	if( .format == "ped" ){
		data.new <- Haplin:::f.ped.to.mfc.new( data.in, design )
	} else {
		data.new <- data.in
		
		if( !is.null( data.new$cov.data ) ){
			orig.colnames <- colnames( data.new$cov.data )
			new.colnames <- paste0( orig.colnames, ".c" )
			colnames( data.new$cov.data ) <- new.colnames
		}
	}
	
	# re-coding the variables per column
	data.recoded <- Haplin:::f.prep.data.parallel( data.new, design, marker.names, ncpu )
	
	## add information
	class( data.recoded ) <- "haplin.ready"
	data.recoded$aux$info <- .info
	data.recoded$aux$class <- class( data.recoded )
	data.recoded$aux$marker.names <- marker.names
	data.recoded$aux$map.filename <- map.file
	
	## find rows with missing data, keep for future reference
	sum.na.list <- lapply( data.recoded$gen.data, function( gen.el.ff ){
		gen.el <- Haplin:::f.extract.ff.numeric( gen.el.ff )
		rowSums( is.na( gen.el ) )
	})
	sum.na <- Reduce( cbind, sum.na.list )
	colnames( sum.na ) <- NULL
	.is.na <- sum.na > 0.1
	rows.with.na <- sum( .is.na )
	data.recoded$aux$nrows.with.missing <- rows.with.na
	data.recoded$aux$which.rows.with.missing <- which( .is.na )
	
	## saving the data in the .RData and .ffData files
	cat( "Saving data... \n" )
	cur.names <- c()
	for( i in 1:length( data.recoded$gen.data ) ){
		cur.name <- paste( get( ".gen.cols.name", envir = Haplin:::.haplinEnv ), i, sep = "." )
		assign( cur.name, data.recoded$gen.data[[i]] )
		cur.names <- c( cur.names, cur.name )
	}
	aux <- data.recoded$aux
	# here, we need to change the format, since now it's 6 columns per 1 marker in the genotype matrix
	aux$info$filespecs$format <- "haplin"
	
	save.list <- c( cur.names, "aux" )
	if( !is.null( data.recoded$cov.data ) ){
		cov.data.in <- data.recoded$cov.data
		save.list <- c( save.list, "cov.data.in" )
	}
	ff::ffsave( list = save.list, file = file.path( dir.out, files.list$file.out.base ) )
	cat( "... saved to files:", files.list$file.out.ff, ", ", files.list$file.out.aux, "\n" )
	return( data.recoded )
}
