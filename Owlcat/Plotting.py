# -*- coding: utf-8 -*-
import matplotlib.pyplot as pyplot
import numpy
import numpy.ma
import math
import sys
import traceback

class PlotCollection (object):
  """PlotCollection plots a collection of data tracks in one plot""";
  def __init__ (self,options):
    self.data = {};
    self.mean = {};
    self.stddev = {};
    self.label = {};
    self.track_counts = {};
    # inter-plot offset -- may be set from outside
    self.offset = 0;

  def num_tracks (self):
    return len(self.data);

  def add_track (self,key,data,mean=None,stddev=None,label=None,count=1):
    self.data[key] = data;
    self.mean[key] = mean if mean is not None else data.mean();
    self.stddev[key] = stddev if stddev is not None else data.std();
    self.label[key] = label if label is not None else str(key);
    self.track_counts[key] = count;

  def get_track_data (self,key):
    return self.data.get(key,None);

  def get_track_count (self,key):
    return self.track_counts.get(key,0);

  # function to make single-page plot
  def make_figure (self,keylist=None,suptitle=None,save=None,
                   dual=True,offset_std=10,
                   figsize=(210,290),dpi=100,papertype='a4',landscape=False):
    # setup sizes
    keylist = keylist or sorted(self.data.keys());
    # compute inter-plot offset, based on the median stddev
    offset = self.offset;
    if not offset:
      stddevs = sorted(self.stddev.itervalues());
      if stddevs:
        offset = self.offset = offset_std*stddevs[len(stddevs)/2];
    # if still not set, force to 1
    if not offset:
      offset = self.offset = 1;
    # partition keys into two sets, if asked to
    if dual:
      # split ifrs into two sets (if number is odd, make sure first subset has the extra 1)
      n = len(keylist);
      n2 = len(keylist)/2;
      if n%2:
        n2 += 1;
      keysets = [ keylist[0:n2],keylist[n2:] ];
    else:
      keysets = [ keylist ];
    # create Figure object of given size and resolution
    figsize_in = (figsize[0]/25.4,figsize[1]/25.4);
    fig = pyplot.figure(figsize=figsize_in,dpi=dpi);
    # margin sizes. We want to keep them fixed in absolute terms,
    # regardless of plot size. The relative numbers here are goof for a 210x290 plot, so
    # we rescale them accordingly.
    mleft   = 0.05 * 210./figsize[0];
    mbottom = 0.01 * 290./figsize[1];
    mright  = 0.01 * 210./figsize[0];
    mtop    = 0.03 * 290./figsize[1];
    ytitle  = 1 - 0.01 * 290./figsize[1];
    width   = ( 1 - mleft - mright );
    height  = ( 1 - mtop - mbottom );
    # now plot all tracks
    for iplot,keys in enumerate(keysets):
      if dual:
        plt = fig.add_subplot(1,2,iplot+1);
      else:
        plt = fig.add_axes([mleft,mbottom,width,height]);
      nx = max([len(self.data.get(key,[])) for key in keys]);
      # first and last key
      key0 = keys[0];
      key1 = keys[-1];
      mindata,maxdata = 0,1;
      # make plots for all ifrs
      firstkey = None;
      for key in keys:
        data = self.data.get(key);
        if data is None or data.mask.all():
          continue;
        else:
          # print data.min(),data.max();
          # first plot goes at its natural coordinates.
          if firstkey is None:
            firstkey = key;
            y0 = 0;
            mindata = y0text = self.mean[key] - offset/2;
            dd = data;
          # subsequent plots offset accordingly
          else:
            y0 += offset;
            dd = data + y0;
          maxdata = self.mean[key] + y0 + offset/2;
          # plot data
          try:
            plt.plot(dd);
          except:
            traceback.print_exc();
            print "Error plotting data for",key;
            continue;
          if self.label[key]:
            # figure out where to place text label
            ytext = dd[:nx/10].mean();
            if numpy.isnan(ytext):
              ytext = dd.mean();
            # do not put too close to previous label
            y0text = max(y0text+offset/4,ytext);
            plt.annotate(self.label[key],xy=(nx/100,ytext),xytext=(nx/100,y0text),
                  size=5,horizontalalignment='left',verticalalignment='center',
                  bbox=dict(facecolor='white',edgecolor='none',alpha=0.6),
                  arrowprops=dict(fc='none',width=0,headwidth=0,ec="0.8",alpha=0.6));
            #plt.text(0,self.mean[key]+y0,self.label[key],
                  #size=5,horizontalalignment='left',verticalalignment='center',
                  #bbox=dict(facecolor='white',edgecolor='none',alpha=0.6));
          
      # set plot limits. First panel is determined by data. Scale of second panel is fixed to first.
      plt.set_xbound(0,nx);
      if iplot == 0:
        ylim = mindata,maxdata;
        plt.set_ybound(*ylim);
        for lab in plt.get_yticklabels():
          lab.set_fontsize(5);
      else:
        plt.set_yticklabels([]);
        plt.set_ybound(*ylim);
      # other plot frills
      plt.set_xticklabels([]);
    fig.subplots_adjust(left=mleft,right=1-mright,top=1-mtop,bottom=mbottom,wspace=0.01);
    # plot title if asked to
    if suptitle:
      fig.suptitle(suptitle,y=ytitle,size=8);
    if save:
      fig.savefig(save,papertype=papertype,
                  orientation='portrait' if not landscape else 'landscape');
      print "===> Wrote",save;
    return fig;

