#' @importFrom basilisk BasiliskEnvironment
env1 <- BasiliskEnvironment("env1", pkgname="dnacycpv2",
                            packages=c('python==3.11',
                                       'pandas==2.1.2',
                                       'tensorflow==2.14.0',
                                       'keras==2.14.0',
                                       'docopt==0.6.2'),
                            pip=c('numpy==1.26.1',
                                  'bio==1.7.1'),
                            path=c('dnacycpv2_python'))

#' @importFrom basilisk BasiliskEnvironment
env2 <- BasiliskEnvironment("env2", pkgname="dnacycpv2",
                            packages=character(0), path="dnacycpv2_python")
