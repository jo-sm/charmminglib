"""
TODO
"""


# Plotting
    def plot(self,timeSeries,title,write=False):
        """
        Plots a figure containing radius of gyration as a function of time.
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        """
        # Plotting
        pyplot.plot(map(lambda x: self.xAxisScaleFactor*x,range(len(timeSeries))),timeSeries,'darkblue')
        pyplot.xlabel(r'$time\ (ns)$')
        pyplot.ylabel(r'$Radius\ of\ Gyration\ (\AA)$')
        pyplot.axis([0,self.xAxisScaleFactor*len(timeSeries),floor(timeSeries.min()),ceil(timeSeries.max())])
        pyplot.title(title)
        # Show/Write
        if write:
            FileName    = '%s%sRgPlot_%s.%s' % (self.anlPath,os.sep,title,write)
            pyplot.savefig(FileName,format=write)
        else:
            pyplot.show()

    def plotTempFig(self,TempList,RgList,title,write=False):
        """
        Plots average Radius of Gyration as a function of temperature.
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        """
        # Plotting
        pyplot.plot(TempList,RgList,'g^')
        pyplot.xlabel(r'$temperature\ (K)$')
        pyplot.ylabel(r'$Radius\ of\ Gyration\ (\AA)$')
        pyplot.title(title)
        # Show/Write
        if write:
            FileName    = '%s%sRgPlot_%s.%s' % (self.anlPath,os.sep,title,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)
        else:
            pyplot.show()
