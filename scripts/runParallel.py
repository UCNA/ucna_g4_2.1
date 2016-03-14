from subprocess import call
from multiprocessing import Process



def runGeant(execname):
    call(['nice', './ucn', macname])


if __name__ == '__main__':

    jobs = []


    ##
    ## Edit the list of executables here
    ## Each entry in the list is a string specifying the executable name
    ## Example is a bunch of execubales called "testi", with 0 <= i <= 9.
    ##
    executables = [ "test0", "test1", "test2", "test3", "test4",
                    "test5", "test6", "test7", "test8", "test9" ]
    
    macfiles = ["run0.mac", "run1.mac", "run2.mac", "run3.mac", "run4.mac"
                "run5.mac", "run6.mac", "run7.mac", "run8.mac", "run9.mac" ] 

    for m in macfiles:
            print 'Running simulation: '+ ' ./ucn ' + m
            jobs.append(Process(target=runGeant, args=('./ucn',)))
            jobs.append(Process(target=runGeant, args=(m,)))

    for j in jobs:
        j.start()

    for j in jobs:
        j.join()

