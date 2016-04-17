# R2WBwrapper.R: a R2WinBUGS wrapper
# KUBO Takuya (kubo@ees.hokudai.ac.jp)
# DATE: 2014-04-25
# This software is distributed under the terms of
# the GNU General Public License Version 2, June 1991.

cat("# reading \"R2WBwrapper.R\" (written by kubo@ees.hokudai.ac.jp)...\n")
library(R2WinBUGS)
library(coda)

setG <- function(name, x) assign(name, x, envir = .GlobalEnv)

# Clear ==================================================================
clear.v.datanames <- function()
{
	setG("v.datanames", NULL)
}
clear.param <- function()
{
	setG("parameters.to.save", NULL)
	setG("list.inits", NULL)
}
clear.data.param <- function()
{
	clear.v.datanames()
	clear.param()
}

# !!! Initialization !!!
clear.data.param()

# Data ===================================================================
set.data.sub <- function(name, val)
{
	if (is.null(val)) stop("# ERROR: no data for", name, "\n")
	ditem <- val
	setG(name, ditem)
}

set.data <- function(name = NA, val = NA, list.data = NA)
{
	if (!is.na(list.data)) {
		v.datanames <- c(v.datanames, names(list.data))
		setG("v.datanames", v.datanames)
		for (name in v.datanames) set.data.sub(name, list.data[[name]])
	} else if (!is.na(name)) {
		v.datanames <- c(v.datanames, name)
		setG("v.datanames", v.datanames)
		if (length(dim(val)) < 2) val <- as.vector(val)
		set.data.sub(name, val)
	} else {
		stop("# Error: set either (name, val)  or list.data")
	}
}

# Parameters =============================================================

# Examples:
# set.param("a", 0, save = FALSE)
# set.param("b", function() rnorm(100, mean = 0, sd = 0.1))
# set.param("c", NA)
set.param <- function(name, val = NA, save = TRUE)
{
	parameters.to.save <- get("parameters.to.save", envir = .GlobalEnv)
	list.inits <- get("list.inits", envir = .GlobalEnv)
	if (save) {
		parameters.to.save <- c(parameters.to.save, name)
	}
	setG("parameters.to.save", parameters.to.save)
	if (any(is.na(val))) return() # do not initialize
	# parameter initialization
	if (class(val) == "numeric" & length(dim(val)) < 2) val <- as.vector(val)
	list.new <- list(val)
	names(list.new) <- name
	list.inits <- c(list.inits, list.new)
	setG("list.inits", list.inits)
}

# Calling WinBUGS ========================================================
bugscalling <- function(
	model.file = "model.bug.txt",
	debug = FALSE,
	n.chains = 3,
	inits,
	n.iter, n.burnin, n.thin,
	...
) {
	WINEPATH <- "/usr/bin/winepath"
	sys.time <- system.time(
		post.bugs <- bugs(
			data = v.datanames,
			inits = inits,
			parameters.to.save = parameters.to.save,
			model.file = model.file,
			n.chains = n.chains, n.iter = n.iter,
			n.burnin = n.burnin, n.thin = n.thin,
			bugs.directory = paste(
				Sys.getenv("HOME"),
				".wine/drive_c/Program\ Files/WinBUGS14/",
				sep = "/"
			),
			#working.directory = NULL,
			clearWD = TRUE,
			useWINE = TRUE,
			newWINE = TRUE,
			WINE = "/usr/bin/wine",
			WINEPATH = WINEPATH,
			debug = debug,
			...
		)
	)
	print(sys.time)
	return(post.bugs)
}

call.bugs <- function(
	file = "model.bug.txt",
	debug = FALSE,
	n.cores = 1, # number of cores, 1: non-parallel
	n.chains = 3,
	n.iter, n.burnin, n.thin,
	bugs.seed = NULL,
	...
) {
	if (any(is.null(bugs.seed))) {
		bugs.seed <- sample(-10^6:10^6, 1)
	}
	model.file <- paste(getwd(), file, sep = "/")
	inits <- function() list.inits
	if (Sys.info()["sysname"] == "Linux") {
		if (n.cores > 1) {
			library(foreach)
			library(doMC)
			registerDoMC(cores = n.cores)
			list.bugs <- foreach(job = 1:n.chains) %dopar% {
				wd <- sprintf("tmp%i", job) # i.e. "tmp1", "tmp2", "tmp3"
				if (!any(dir() %in% wd)) dir.create(wd)
				bugscalling(
					model.file = model.file,
					n.chains = 1,
					debug = FALSE,
					inits = inits,
					bugs.seed = bugs.seed,
					working.directory = wd,
					n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin
				)
			}
			post.bugs <- merge.bugs(list.bugs)
		} else {
			post.bugs <- bugscalling(
				model.file = model.file, n.chains = n.chains, debug = debug,
				inits = inits,
				n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin,
				...
			)
		}
	} else {
		post.bugs <- bugs(
			v.datanames, inits, parameters.to.save, model.file,
			n.chains = n.chains, n.iter = n.iter,
			n.burnin = n.burnin, n.thin = n.thin,
			bugs.directory = "c:/Program Files/WinBUGS14/",
			working.directory = NULL,
			clearWD = TRUE,
			debug = debug,
			bugs.seed = bugs.seed,
			...
		)
	}
	clear.data.param() # !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!!
	post.bugs
}

# Post processings =======================================================
to.list <- function(x.bugs) mcmc.list(
	lapply(
		1:x.bugs$n.chain,
		function(chain) as.mcmc(x.bugs$sims.array[, chain,])
	)
)

to.mcmc <- function(x.bugs) as.mcmc(x.bugs$sims.matrix)

pl <- function(key, x.bugs = post.bugs, ...)
{
	cn <- colnames(to.mcmc(x.bugs))
	if (!exists("post.list")) post.list <- to.list(x.bugs)
	if (length(key) == 1) {
		plot(post.list[,grep(key, cn),], ask = T, smooth = F, ...)
	} else {
		plot(post.list[,cn %in% key,], ask = T, smooth = F, ...)
	}
}

pg <- function(x.bugs = post.bugs, digits.summary = 3, ...) page(
	x.bugs, "print", digits.summary = digits.summary, ...
)

choose.post <- function(key, alpha = 0.05, x.bugs = post.bugs)
{
	x.mcmc <- to.mcmc(x.bugs)
	c.mcmc <- x.mcmc[, grep(key, colnames(x.mcmc))]
	for (i in colnames(c.mcmc)) {
		v.prob <- quantile(
			c.mcmc[, i],
			probs = c(alpha * 0.05, 0.50, 1 - alpha * 0.5)
		)
		if (v.prob[1] * v.prob[3] > 0) {
			cat(sprintf(
				"%s\t%8.2f%8.2f%8.2f\n",
				i, v.prob[1], v.prob[2], v.prob[3]
			))
		}
	}
}

# car.normal() ===========================================================
get.cn.parameters <- function(vx, vy, vid)
{ # min(vx) == 1 and min(vy) == 1
	m <- matrix(0, max(vx) + 2, max(vy) + 2)
	n <- length(vid)
	for (i in 1:n) m[vx[i] + 1, vy[i] + 1] <- vid[i]
	v.range <- c(-1, 0, 1) + 1 # 8 neighbors + self (0, 0)
	neighbor.id <- sapply(
		1:n, function(i) c(m[vx[i] + v.range, vy[i] + v.range])
	)[-5,] # to remove self
	s.zero <- c(neighbor.id) == 0
	list(
		Adj = c(neighbor.id)[!s.zero],
		Num = apply(neighbor.id, 2, function(X) sum(X > 0))
	)
}

# merge.bugs() for parallel run using libray(foreach)
merge.bugs <- function(list.bugs)
{
	b1 <- list.bugs[[1]]
	m1 <- b1$sims.matrix
	mall <- matrix(numeric(0), 0, ncol(m1))
	n.chains <- length(list.bugs)
	for (i in 1:n.chains) {
		mall <- rbind(mall, list.bugs[[i]]$sims.matrix)
	}
	sims.array <- array(mall, dim = c(nrow(m1), n.chains, ncol(m1)))
	dimnames(sims.array) <- list(NULL, NULL, colnames(m1))
	as.bugs.array(
		sims.array = sims.array,
		model.file = b1$model.file,
		program = b1$program,
		DIC = TRUE,
		DICOutput = NULL,
		n.iter = b1$n.iter,
		n.burnin = b1$n.burnin,
		n.thin = b1$n.thin
	)
}

