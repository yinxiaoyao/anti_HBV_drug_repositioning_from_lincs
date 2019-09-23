# helper function to subset down to landmark probes only
mkgct = function(ofile,ds,precision=4,appenddim=T,ver=3) {
  # gct must contain the following fields
  #          mat: Numeric data matrix [RxC]
  #          rid: Cell array of row ids
  #          rhd: Cell array of row annotation fieldnames
  #          rdesc: Cell array of row annotations
  #          cid: Cell array of column ids
  #          chd: Cell array of column annotation fieldnames
  #          cdesc: Cell array of column annotations
  #          version: GCT version string
  #          src: Source filename
  
  
  # append the dimensions of the data set, if desired
  nc = ncol(ds@mat)
  nr = nrow(ds@mat)
  if (appenddim) {
    outFile = basename(ofile)
    filename = strsplit(outFile,'.',fixed=T)[[1]][1]
    ofile = path.join(dirname(ofile),
                      sprintf('%s_n%dx%d.gct',filename,
                              nc,nr)) 
  }
  
  precision = floor(precision)
  cat(sprintf('Saving file to %s\n',ofile))
  cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
  cat(sprintf('Setting precision to %d\n',precision))
  
  # open file      
  if (ver==3) {
    nrdesc = dim(ds@rdesc)[2]
    ncdesc = dim(ds@cdesc)[2]
    colkeys = colnames(ds@cdesc)
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id',colnames(ds@rdesc),ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    # line 4 + ncdesc: sample desc
    filler = 'na'
    for (ii in 1:ncdesc) {
      if (is.numeric(ds@cdesc[,ii])) {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    round(ds@cdesc[,ii],precision)),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)  
      } else {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    ds@cdesc[,ii]),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)
      }
    }
  } else {
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id','Description',ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
  }
  
  for (ii in 1:nr) {    
    # print rows
    cat(paste(c(ds@rid[ii],
                ds@rdesc[ii,],
                round(ds@mat[ii,],precision)),collapse='\t'),
        sep='\n',file=ofile,append=T)
  }
  cat(sprintf('Saved.\n'))  
}

subset.to.landmarks <- function(ds) {
 options(stringsAsFactors=FALSE)
 lm <- read.table('/xchip/cogs/data/vdb/spaces/lm_epsilon_n978.grp') 
 return(gctextract.tool(ds,rid=lm$V1))
}

# helper function to set all the row and column annotations to the correct data type
fix.datatypes = function(meta) {
    # turn all warnings to errors so we can use the try statement to grab strings
    options(warn = 2)
    for (field.name in names(meta)) {
        # get the field
        field = meta[[field.name]]
        # check if it's numeric
        try({field = as.numeric(field)}, silent = TRUE)
        if (is.numeric(field)) {
            #check if it's an integer
            int.field = NULL
            try({int.field = as.integer(field)}, silent = TRUE)
            if ( ! is.null(int.field) && identical(int.field, field) )
                field = int.field
        }
        # insert back into the annotations
        meta[[field.name]] = field
    }
    options(warn = 0)
    return(meta)
}

# define the gct object class
setClass("GCT",
         representation(
             mat = "matrix",
             rid = "vector",
             cid = "vector",
             rdesc = "data.frame",
             cdesc = "data.frame",
             version = "character",
             src = "character"
         )
)

# define the initialization method for the class
setMethod("initialize",
          signature = "GCT",
          definition = function(.Object, src, rid = NULL, cid = NULL) {
              # check to make sure it's either .gct or .gctx
              if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                  stop("Either a .gct or .gctx file must be given")
              if (grepl(".gct$", src)) {
                  if ( ! is.null(rid) || !is.null(cid) )
                      stop("rid and cid values may only be given for .gctx files, not .gct files")
                  # parse the .gct
                  .Object@src = src
                  # get the .gct version by reading first line
                  .Object@version = scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
                  # get matrix dimensions by reading second line
                  dimensions = scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
                  nrmat = dimensions[1]
                  ncmat = dimensions[2]
                  nrhd = dimensions[3]
                  nchd = dimensions[4]
                  # read in header line
                  header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
                  # construct row header and column id's from the header line
                  if ( nrhd ) {
                      rhd = header[2:(nrhd+1)]
                      cid = header[-(nrhd+1):-1]
                  }
                  else {
                      rhd = NULL
                      cid = header[-1]
                  }
                  # read in the next set of headers (column annotations) and shape into a matrix
                  if ( nchd ) {
                      header = scan(src, what = "", nlines = nchd, skip = 3, sep = "\t", 
                                    quote = NULL, quiet = TRUE)		
                      header = matrix(header, nrow = nchd, 
                                      ncol = ncmat + nrhd + 1, byrow = TRUE)
                      # extract the column header and column descriptions
                      chd = header[,1]
                      cdesc = header[,-(nrhd+1):-1]
                      # need to transpose in the case where there's only one column annotation
                      if ( nchd == 1 )
                          cdesc = t(cdesc)
                  }
                  else {
                      chd = NULL
                      cdesc = data.frame()
                  }
                  # read in the data matrix and row descriptions, shape into a matrix
                  mat = scan(src, what = "", nlines = nrmat, 
                             skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
                  mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + 1, 
                               byrow = TRUE)
                  # Extract the row id's row descriptions, and the data matrix
                  rid = mat[,1]
                  if ( nrhd ) {
                      # need as.matrix for the case where there's only one row annotation
                      rdesc = as.matrix(mat[,2:(nrhd + 1)])
                      mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                                   nrow = nrmat, ncol = ncmat)
                  }
                  else {
                      rdesc = data.frame()
                      mat = matrix(as.numeric(mat[,-1]),
                                   nrow = nrmat, ncol = ncmat)
                  }
                  # assign names to the data matrix and the row and column descriptions
                  dimnames(mat) = list(rid, cid)
                  if ( nrhd ) {
                      dimnames(rdesc) = list(rid,rhd)
                      rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
                  }
                  if ( nchd ) {
                      cdesc = t(cdesc)
                      dimnames(cdesc) = list(cid,chd)
                      cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
                  }
                  # assign to the GCT slots
                  .Object@mat = mat
                  .Object@rid = rownames(mat)
                  .Object@cid = colnames(mat)
                  .Object@rdesc = fix.datatypes(rdesc)
                  .Object@cdesc = fix.datatypes(cdesc)
                  return(.Object)
              }
              else {
                  # parse the .gctx
                  library(h5r)
                  .Object@src = src
                  # if the rid's or column id's are .grp files, read them in
                  if ( length(rid) == 1 && grepl(".grp$", rid) )
                      rid = parse.grp(rid)
                  if ( length(cid) == 1 && grepl(".grp$", cid) )
                      cid = parse.grp(cid)
                  # connect to the file; get row and column id's if needed
                  f = H5File(src)
                  .Object@version = getH5Attribute(f, "version")[]
                  mat.data = getH5Dataset(getH5Group(f, "0/DATA/0"), "matrix")
                  row.meta = getH5Group(f, "0/META/ROW")
                  col.meta = getH5Group(f, "0/META/COL")
                  # get the row id's and cid's from the file
                  rid.data = getH5Dataset(row.meta, "id")
                  rid.all = rid.data[]
                  cid.data = getH5Dataset(col.meta, "id")
                  cid.all = cid.data[]
                  # if rid's and cid's aren't given by user, assign them
                  if ( is.null(rid) )
                      rid = rid.all
                  if ( is.null(cid) )
                      cid = cid.all
                  # get the row and column indices we need
                  ridx = match(rid, rid.all)
                  if ( any(is.na(ridx)) ) {
                      missing.rid = sort(rid[is.na(ridx)])
                      missing.str = paste(missing.rid, collapse = "\n")
                      warning(sprintf("There were no matches for the following rid's:\n%s", missing.str))
                      rid = rid[!is.na(ridx)]
                      ridx = ridx[!is.na(ridx)]
                  }
                  cidx = match(cid, cid.all)
                  if ( any(is.na(cidx)) ) {
                      missing.cid = sort(cid[is.na(cidx)])
                      missing.str = paste(missing.cid, collapse = "\n")
                      warning(sprintf("There were no matches for the following cid's:\n%s", missing.str))
                      cid = cid[!is.na(cidx)]
                      cidx = cidx[!is.na(cidx)]
                  }
                  # read the data matrix
                  if ( length(cidx) == 1)
                      mat = as.matrix(mat.data[cidx, ridx])
                  else
                      mat = t(mat.data[cidx, ridx])
                  dimnames(mat) = list(rid, cid)
                  
                  # read the row meta data
                  row.contents = listH5Contents(row.meta)
                  row.fields = names(row.contents)
                  row.fields = row.fields[row.fields != "." & row.fields != "id"]
                  # preallocate the data frame
                  rdesc = data.frame(matrix(rep("", length(ridx) * length(row.fields)),
                                            nrow = length(ridx), ncol = length(row.fields),
                                            dimnames = list(rid, row.fields)))
                  # read it in
                  for ( row.field in row.fields ) {
                      row.data = getH5Dataset(row.meta, row.field)
                      rdesc[[row.field]] = row.data[ridx]
                  }
                  
                  # read the column meta data
                  col.contents = listH5Contents(col.meta)
                  col.fields = names(col.contents)
                  col.fields = col.fields[col.fields != "." & col.fields != "id"]
                  # preallocate
                  cdesc = data.frame(matrix(rep("", length(cidx) * length(col.fields)),
                                            nrow = length(cidx), ncol = length(col.fields),
                                            dimnames = list(cid, col.fields)))
                  # read in
                  for ( col.field in col.fields ) {
                      col.data = getH5Dataset(col.meta, col.field)
                      cdesc[[col.field]] = col.data[cidx]
                  }
                  
                  # assign all data
                  .Object@mat = mat
                  .Object@rid = rid
                  .Object@cid = cid
                  .Object@rdesc = fix.datatypes(rdesc)
                  .Object@cdesc = fix.datatypes(cdesc)
                  return(.Object)
              }
          }
)

# function to parse a GCT(X)
# just instantiates a new .gct object
parse.gctx = function(fname, rid = NULL, cid = NULL) {
    ds = new("GCT", src = fname, rid = rid, cid = cid)
    return(ds)
}

# method for gct class object to extract rows and cols from .gct file
gctextract.tool = function(ds, rid = NULL, cid = NULL) {
    if (! is.null(rid)) {
        # these will be ordered as rid (that's how the intersect works in R)
        rid.ds = intersect(rid, ds@rid)
        if (! identical(rid.ds, rid)) {
            missings = setdiff(rid, rid.ds)
            warning("The following rid's were not found: ", 
                    paste(missings, collapse = ", "))
        }
        ds@mat = ds@mat[rid.ds,,drop = FALSE]
        ds@rid = rid.ds
        ds@rdesc = ds@rdesc[rid.ds,]
    }
    # extract the columns
    if (! is.null(cid)) {
        # will be ordered as cid
        cid.ds = intersect(cid, ds@cid)
        if (! identical(cid.ds, cid)) {
            missings = setdiff(cid, cid.ds)
            warning("The following cid's were not found: ",
                    paste(missings, collapse = ", "))
        }
        ds@mat = ds@mat[,cid.ds,drop = FALSE]
        ds@cid = cid.ds
        ds@cdesc = ds@cdesc[cid.ds,]
    }
    return(ds)
}

# function to read .grp files
parse.grp = function(fname) {
    grp = scan(fname, what = "", quote = NULL, quiet = TRUE)
    return(grp)
}

# function to read .gmx files
parse.gmx = function(fname) {
    tmp = read.table(fname, sep = "\t", 
                     header = TRUE, stringsAsFactors = FALSE)
    # preallocate a list for the gmx
    L = list()
    # loop over the first row of the .gmx
    for ( n in names(tmp) ) {
        # get all the values; remove empties at the end
        values = tmp[[n]][-1]
        remove.idx = values == ""
        values = values[!remove.idx]
        # put in a list
        L[[n]] = list(head = n,
                      desc = tmp[[n]][1], 
                      len = length(values), 
                      entry = values)
    }
    return(L)
}

# function to read .gmt files
parse.gmt = function(fname) {
    gmt.lines = scan(fname, what = "", sep = "\n",
                     quote = NULL, quiet = TRUE)
    tmp = lapply(gmt.lines, function(x) unlist(strsplit(x, "\t")))
    mk.gmt.entry = function(x) {
        L = list()
        L[["head"]] = x[1]
        L[["desc"]] = x[2]
        l.entry = x[-c(1:2)]
        idx = l.entry != ""
        L[["entry"]] = l.entry[idx]
        L[["len"]] = length(L[["entry"]])
        return(L)
    }
    L = lapply(tmp, function(x) mk.gmt.entry(x))
    names(L) = unlist(lapply(L, function(x) x$head))
    return(L)
}

# function to write tab-delimited text files at a fixed numerical precision
mktbl = function(tbl, ofile, precision = 4, col.names = TRUE, row.names = TRUE) {
    tbl = as.data.frame(tbl)
    # format numeric floats to two decimal points; leave all others as is
    for (col in names(tbl))
        if (class(tbl[[col]]) == "numeric")
            if (! (all(round(tbl[[col]]) == tbl[[col]]))) {
                # the format string; tells us the precision
                fmt = paste("%0.", as.character(precision), "f", sep = "")
                tbl[[col]] = sprintf(fmt, tbl[[col]])
            }
    write.table(tbl, file = ofile, sep = "\t", quote = FALSE, 
                col.names = col.names, row.names = row.names)
}

# function to join a bunch of arguments into a file path (same as matlab's fullfile)
path.join = function(...) {
    args = c(...)
    args = sub("/$", "", args)
    path = paste(args, collapse = "/")
    return(path)
}

# function to make a RESTful call to the supplied URL
# and return the parsed JSON response
call.api <- function(url) {
  require(RCurl)
  require(rjson)
  # call API, parse JSON, and return parsed object
  # .opts = list(ssl.verifypeer = FALSE)
  # disables the SSL certificate check
  return(fromJSON(getURL(url, .opts = list(ssl.verifypeer = FALSE))))
}