'''
# ------------------------------------------------------------------------------
#
# Created by : Mary Clay PhD
# e-mail: mary.clay@stjude.org
# St Jude Children's Research Hospital 
# Department of Structural Biology Memphis, TN 
# 11/10/2022
#
# Updated January 30, 2023 to include peak intensities 
# Updated March 09, 2023 to account for old cycle7.peaks format missing information for multiple assignment crosspeaks
# Updated March 14, 2023 to add header to output list, and short constraints tag
# Updated July 17, 2023 to include heavy atom distances to assess possible diffusion peaks
# ------------------------------------------------------------------------------

'''
import os
import sys
import numpy as np
import re
import glob

AAA_dict = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "HIST": "H","HIS+": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": 'V', "MSE":'M', "PTR":'Y', "TPO":"T", "SEP":'S'}


def analize_noa(Seqdict, violdict, qupldict, upldict, pad, upldict2, distDF,fillteredDF, posdiffout):
  cwd = os.getcwd() + '/'
  calc = cwd + 'CALC.cya'
  outdir = cwd +'noa_analysis/'
  noa7 = cwd +'cycle7.noa'
  cya_plists = [line.strip() for line in open(calc).readlines() if line.strip() and 'peaks' in line and not re.match(r'^\s*#', line)][0].split()[2].split(',')
  prots = [line.strip() for line in open(calc).readlines() if line.strip() and 'prot' in line and not re.match(r'^\s*#', line)][0].split()[2].split(',')
  log = glob.glob(os.path.join(cwd + 'log*'))[0]
  header = [ '### UPL: Peak was idenitified in the output.upl file\n',
  '### swapped: steriospecific assignment was swapped relative to the input chemical shift list\n'
  '### unsed: Peaks which CYANA did not use assignment possiblity #unused\n',
  '### no assignment: Peaks which CYANA found no assignment possibility #no assignmnet\n',
  '### low suport: Peaks with SUP < 0.6 \n',
  '### poor chem shit: Peaks with poor chemical shift matching (Pshift < 0.6)\n',
  '### long distance: Peak associated with UPL > 6.0\n',
  '### short distance: Peak associated with UPL < 3.0 \n',
  '### lone: Peak represening the only instance of an assignment\n',
  "#{:}  {:^26}  {:^24}  {:^5}  {:^10}  {:^6}  {:^24}\n".format('Peak #','Frequencies','Connection','UPL', 'Range' ,'Pshift','Comment')
  ]
  plist_dict = {}
  for x in range(len(cya_plists)):
    plist = cya_plists[x].replace('.peaks','')
    plist_dict[plist] = str(x)
    exec("outlist{:} = []".format(str(x)))
    exec("usedquestlist{:} = []".format(str(x)))
    exec("peaks{:} = {{}}".format(str(x)))
    exec("intensity{:} = {{}}".format(str(x)))
  for x in range(len(cya_plists)):
    intdict = eval('intensity' + str(x))
    pdict = eval('peaks' + str(x))
    peaklines = open(cwd + cya_plists[x].replace('.peaks','-cycle7.peaks')).readlines()
    ndim = int(peaklines[0][-2])
    intidx = ndim + 3
    VCV = 2*ndim + 8
    for i in range(len(peaklines)):
      line = peaklines[i]
      if line.strip():
        if not re.match(r'^\s*#', line):
          if line[0:7] == '       ':  # account for old format of cycle7.peaks
            line = peaklines[i-1][0:peaklines[i-1].find(' U ')+35] + ' ' + peaklines[i].strip()
            peaklines[i] = line
          pdict[int(line.split()[0])] = [float(line.split()[1]),float(line.split()[2]),float(line.split()[3])]
          if re.search('#VC',line):
            VC = float(line.split()[VCV])
          else: VC = 1.0
          intensity = "{:12.6E}".format(float(line.split()[intidx]) *VC)
          if int(line.split()[0]) in intdict.keys():
            intdict[int(line.split()[0])].append(intensity)
          else:
            intdict[int(line.split()[0])] = [intensity]

  Calibration_cns = []
  Swapped = {}
  for line in open(log).readlines():
    if "Calibration constant for peak list" in line:
      newline = line.replace(line.split()[5],cya_plists[int(line.split()[5][0:-1]) -1].replace('.peaks',':'))
      if newline not in Calibration_cns:
        Calibration_cns.append(newline)
    if line.strip() and line.strip().split()[-1] == 'swapped':
      res = line.strip().split()
      group1 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[2])
      group2 = '{:}{:}-{:}'.format(AAA_dict[res[1]],res[0], res[3])
      Swapped[group1] = group2
      Swapped[group2] = group1
  Assignments= []
  for prot in prots:
    for line in open(cwd + prot.replace('.prot','-final.prot')).readlines():
      if line.strip():
        if '#' != line.strip()[0]:
          resi = line.strip().split()[4]
          atom = line.strip().split()[3]
          assign = Seqdict[resi] + '-' + atom
          if assign not in Assignments and atom[0] in ['H','Q']:
            Assignments.append(assign)
  Assignments = sorted(Assignments, key = lambda x: (x.split('-')[0][1:], x.split('-')[1]))

  assigndict,assigndict2,notassigndict = {}, {},{}
  from itertools import combinations
  ADpairs = [ '{:}-{:}'.format(comb[0],comb[1]) for comb in combinations(Assignments,2)]
  for comb in ADpairs:
    assigndict[comb] = []
    assigndict2[comb] = []
    notassigndict[comb] = []
  noalines = open(noa7).readlines()
  for x in range(len(noalines)):
    line = noalines[x]
    if 'Peak' in noalines[x] and 'out of' in noalines[x+1]:
      if '0 out of' not in noalines[x+1] and 'diagonal' not in noalines[x]:
        peak = int(noalines[x].strip().split()[1])
        plist = noalines[x].strip().split()[3].replace('.peaks','')
        pdict = eval('peaks' + plist_dict[plist])
        intdict = eval('intensity' + plist_dict[plist])
        outlist = eval('outlist' + plist_dict[plist])
        used =eval('usedquestlist' + plist_dict[plist])
        SUP = noalines[x+1].split()[-1].replace(':','')
        linepad = pad[len(plist):]
        Calconst = float(Calibration_cns[int(plist_dict[plist])].split()[-1])
        # print(plist)
        # print('{:1.3E}'.format(Calconst))
        for y in range(2,int(noalines[x+1].split()[0])+2,1):
          cns = noalines[x+y].strip().split()
          if not re.match(r'^\s+[A-Z]+',noalines[x+y]):
            atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5], cns[6], int(cns[7]), float(cns[10])/100, cns[13]
          if re.match(r'^\s+[A-Z]+',noalines[x+y]):
            atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[0],cns[1],int(cns[2]), cns[4], cns[5], int(cns[6]), float(cns[9])/100, cns[12]
          group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
          group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
          note = ''
          if group1 in Swapped.keys(): group1 = Swapped[group1]; note = note + 'swapped '
          if group2 in Swapped.keys(): group2 = Swapped[group2]; note = note + 'swapped '
          conect = '{:}-{:}'.format(group1,group2)
          if '{:} peak {:} from {:}'.format(conect,peak,plist) in upldict2:
            note = 'UPL ' + note
          if float(intdict[peak][y-2]) != 0.0: 
            dist = (Calconst/float(intdict[peak][y-2]))**(1/6)
            drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
          if float(intdict[peak][y-2]) == 0.0: 
            print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
          outline = '{:^28} {:^14} {:>9}A  Peak {:4} from {:<}{:}  pshift {:3.2f} {:}\n'.format(conect,intdict[peak][y-2],drange,peak,linepad,plist,pshift,note)
          udist = ' '
          if conect in upldict.keys(): udist = upldict[conect] + 'A'
          if '{:}-{:}'.format(group1,group2)in ADpairs:
            assigndict2['{:}-{:}'.format(group1,group2)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
            assigndict['{:}-{:}'.format(group1,group2)].append(outline)
          if '{:}-{:}'.format(group2,group1)in ADpairs:
            assigndict['{:}-{:}'.format(group2,group1)].append(outline)
            assigndict2['{:}-{:}'.format(group2,group1)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
          if float(SUP) <= 0.6:
            note = note + 'low support '
            if '{:} {:}'.format(peak, conect) not in used: used.append('{:} {:}'.format(peak, conect))
          if conect in violdict.keys():
            note = note + violdict[conect].replace(' #','').replace("\n",'') + ' '
            if '{:} {:}'.format(peak, conect) not in used: used.append('{:} {:}'.format(peak, conect))
          if conect in qupldict.keys():
            note = note + qupldict[conect]
            if '{:} {:}'.format(peak, conect) not in used: used.append('{:} {:}'.format(peak, conect))
          if pshift <= 0.60 and '{:} {:}'.format(peak, conect) not in used:
            note = note + 'poor chem shift '
            if '{:} {:}'.format(peak, conect) not in used: used.append('{:} {:}'.format(peak, conect))
          outlist.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,note))
          if '{:} {:}'.format(peak, conect) not in used: used.append('{:} {:}'.format(peak, conect))

  print('finished assigned')
  pdiffcount,unassignedcount = 0,0
  ADpairs2 = []
  UNAssigned = []
  for con in ADpairs:
    if len(assigndict[con]) >= 1:
      ADpairs2.append(con)
    if len(assigndict[con]) == 0:
      UNAssigned.append(con)
  for x in range(len(noalines)):
    line = noalines[x]
    if 'Peak' in noalines[x] and '0 out of' in noalines[x+1]:
      peak = int(noalines[x].strip().split()[1])
      plist = noalines[x].strip().split()[3].replace('.peaks','')
      pdict = eval('peaks' + plist_dict[plist])
      used =eval('usedquestlist' + plist_dict[plist])
      intdict = eval('intensity' + plist_dict[plist])
      outlist = eval('outlist' + plist_dict[plist])
      linepad = pad[len(plist):]
      Calconst = float(Calibration_cns[int(plist_dict[plist])].split()[-1])
      if '0 out of 0' not in noalines[x+1]:
        unassignedcount+=1
        localpdif = 0
        probdiffoutlines = []
        nopt = int(noalines[x+1].split()[3])
        for y in range(2,int(noalines[x+1].split()[3])+2,1):
          cns = noalines[x+y].strip().split()
          # if len(cns) > 8:
          if not re.match(r'^\s+[A-Z]+',noalines[x+y]):
            atom1,resn1, resi1, atom2, resn2, resi2, pshift, drange = cns[1],cns[2],int(cns[3]), cns[5],cns[6],int(cns[7]), float(cns[10])/100, cns[13]
          if re.match(r'^\s+[A-Z]+',noalines[x+y]):
            atom1,resn1, resi1, atom2, resn2, resi2,pshift, drange = cns[0],cns[1],int(cns[2]), cns[4],cns[5],int(cns[6]), float(cns[9])/100, cns[12]
          group1 = '{:}{:}-{:}'.format(AAA_dict[resn1],resi1, atom1)
          group2 = '{:}{:}-{:}'.format(AAA_dict[resn2],resi2, atom2)
          if group1 in Swapped.keys(): group1 = Swapped[group1]
          if group2 in Swapped.keys(): group2 = Swapped[group2]
          conect = '{:}-{:}'.format(group1,group2)
          try:
            d = float(distDF.loc[group1,group2].split()[0])
            note = 'unused {:} '.format(noalines[x+1].split()[3])
            common = fillteredDF.dropna(subset=[group1,group2]).index.tolist()
            if d >= 12.0:
              pass
            else: 
              if d >= 8.0 and len(common) > 0:
                common.extend([group1,group2])
                probdiffoutlines.append('{:}\nlong distance heavy {:}\n{:}\n'.format(conect,distDF.loc[group1,group2],fillteredDF.loc[common,common].to_string()))
                # print(conect)
                # print(distDF.loc[group1,group2])
                # print(fillteredDF.loc[common,common])
                note = note + 'prob diff {:}A '.format(distDF.loc[group1,group2])
                localpdif+=1
              if float(intdict[peak][0]) != 0.0: 
                dist = (Calconst/float(intdict[peak][0]))**(1/6)
                drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
              if float(intdict[peak][0]) == 0.0: 
                print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
              outline = '#{:^28} {:^14} {:>9}A  Peak {:4} from {:<}{:}  pshift {:0.2f} unused\n'.format(conect,intdict[peak][0],drange,peak,plist,linepad,pshift)
              if '{:}-{:}'.format(group1,group2)in ADpairs2:
                assigndict['{:}-{:}'.format(group1,group2)].append(outline)
                assigndict2['{:}-{:}'.format(group1,group2)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
              if '{:}-{:}'.format(group2,group1)in ADpairs2:
                assigndict['{:}-{:}'.format(group2,group1)].append(outline)
                assigndict2['{:}-{:}'.format(group2,group1)].append('{:^28}  Peak {:4} from {:<}{:}\n'.format(conect,peak,linepad,plist))
              outline = '#{:^28} {:^14} {:>9}A  Peak {:4} from {:<}{:}  pshift {:0.2f} {:} out of {:} unused {:}\n'.format(conect,intdict[peak][0],drange,peak,plist,linepad,pshift,y-1,noalines[x+1].split()[3],noalines[x+1].split()[3])
              if '{:}-{:}'.format(group1,group2)in UNAssigned:
                notassigndict['{:}-{:}'.format(group1,group2)].append(outline)
              if '{:}-{:}'.format(group2,group1)in UNAssigned:
                notassigndict['{:}-{:}'.format(group2,group1)].append(outline)
              if conect in upldict.keys(): udist = upldict[conect]+'A'
              if conect not in upldict.keys(): udist =' '
              outlist.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6.2f}   {:}\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2],conect,udist, drange +'A',pshift,note))
              used.append('{:} {:}'.format(peak, conect))
          except AttributeError:
            pass
        if localpdif >=1:
          pdiffcount+=1
          posdiffout.write(noalines[x])
          posdiffout.writelines(probdiffoutlines)
      if '0 out of 0' in noalines[x+1]:
        drange = ''
        if float(intdict[peak][0]) != 0.0: 
          dist = (Calconst/float(intdict[peak][0]))**(1/6)
          drange = '{:3.2f}-{:3.2f}'.format(dist, dist*1.25)
        if float(intdict[peak][0]) == 0.0: 
          print('Warning Peak {:>4} from {:} has zero intensity !'.format(peak, plist))
        pdict = eval('peaks' + plist_dict[plist])
        outlist.append("{:>6}  {:>8.3f} {:>8.3f} {:>8.3f}  {:^24}  {:^5}  {:^10}  {:^6}   no assignmnet\n".format(peak,pdict[peak][0],pdict[peak][1],pdict[peak][2], 'none','', drange+'A','na'))
  print('{:} probable diffusion peaks identified in {:} unused peaks'.format(pdiffcount,unassignedcount))
  assigned = open(outdir + 'Assignment_Summary.txt','w')
  unused = open(outdir + 'Unused_Connections.txt','w')
  inconcount,numcon, intracount, shortcount,mediumcount,longcount = 0, 0, 0, 0, 0, 0
  Gamidecount, Tamidecount,Mamidecount = 0, 0,0
  for val in Calibration_cns:
    assigned.write(val.strip()+'\n')
  assigned.write('\n\n')
  for con in ADpairs:
    if len(assigndict[con]) >= 1:
      numcon += 1
      dist = []
      for line in assigndict[con]:
        if 'unused' not in line: dist.append(float(line.split()[2].split('-')[0]))
        if 'unused' in line: dist.append(float(line.split()[3].split('-')[0]))
      if np.std(dist) > 0.3: 
        inconcount += 1
        if con in upldict.keys():
          dref = float(upldict[con])
          if dref < np.min(dist): dref = np.min(dist)
          assigned.write('{:}  {:3.2f}A ({:}): {:0.2f} inconsistent distance\n'.format(con,float(upldict[con]),len(assigndict[con]),np.std(dist)))
        if con not in upldict.keys():
          dref = np.mean(dist)
          assigned.write('{:} ({:}): {:0.2f} inconsistent distance\n'.format(con,len(assigndict[con]),np.std(dist)))
        if np.mean(dist) > dref: ## there are more longer distaces in the set than shorter
          for pline,d in zip(assigndict[con],dist):
            if d <= dref:
              questout =  eval('outlist' + plist_dict[pline.replace("#",'').split()[6]])
              used =eval('usedquestlist' + plist_dict[pline.replace("#",'').split()[6]])
              index = used.index('{:} {:}'.format(int(pline.replace("#",'').split()[4]), pline.replace("#",'').split()[0]))
              pline2 = questout[index]
              questout[index] = pline2.replace('\n', 'inconsistent \n')
              assigned.write(pline.replace('\n',' *\n'))
            else:
              assigned.write(pline)
        assigned.write('\n')
        if np.mean(dist) <= dref:
          for pline,d in zip(assigndict[con],dist):
            if abs(d - np.mean(dist)) > np.std(dist):
              questout =  eval('outlist' + plist_dict[pline.replace("#",'').split()[6]])
              used =eval('usedquestlist' + plist_dict[pline.replace("#",'').split()[6]])
              index = used.index('{:} {:}'.format(int(pline.replace("#",'').split()[4]), pline.replace("#",'').split()[0]))
              pline2 = questout[index]
              questout[index] = pline2.replace('\n', ' inconsistent\n')
              assigned.write(pline.replace('\n','*\n'))
            else:
              assigned.write(pline)
        assigned.write('\n')
      if np.std(dist) < 0.3: 
        if con in upldict.keys():
          assigned.write('{:}  {:3.2f}A ({:}):\n'.format(con,float(upldict[con]),len(assigndict[con])))
        if con not in upldict.keys():
          assigned.write('{:} ({:}):\n'.format(con,len(assigndict[con])))
        assigned.writelines(assigndict[con])
        assigned.write('\n')
    if len(notassigndict[con]) >= 1:
      unused.write('{:} ({:}):\n'.format(con,len(notassigndict[con])))
      unused.writelines(notassigndict[con])
      unused.write('\n')
      rn1 = con.split('-')[0][0]
      rn2 = con.split('-')[2][0]
      if abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) == 0: intracount+= 1
      if abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) <= 4 and abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) != 0: shortcount+= 1
      if abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) <= 5 and abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) >= 2: mediumcount+= 1
      if abs(float(con.split('-')[0][1:]) - float(con.split('-')[2][1:])) > 5 : longcount+= 1
    if len(assigndict[con]) == 1:
      pline = assigndict[con][0].split()
      peak = int(pline[4])
      plist = pline[6]
      pshift = pline[8]
      conect = pline[0]
      pdict = eval('peaks' + plist_dict[plist])
      questout =  eval('outlist' + plist_dict[plist])
      used =eval('usedquestlist' + plist_dict[plist])
      index = used.index('{:} {:}'.format(peak, conect))
      pline = questout[index]
      questout[index] = pline.replace('\n', 'lone\n')


  for x in range(len(cya_plists)):
    exec("good{:} = open('{:}{:}.list','w')".format(str(x), outdir, cya_plists[x].replace('.peaks','')))
    exec("good{:}.writelines(header)".format(str(x)))
    eval("outlist{:}".format(str(x)))
    outlist = eval('outlist' + plist_dict[cya_plists[x].replace('.peaks','')])
    olist = sorted(outlist, key = lambda x: float(x.strip().split()[0]))
    eval("good{:}.writelines(olist)".format(str(x)))
    eval("good{:}.close()".format(str(x)))
  assigned.close()
  # print('Fount {:} inconsistant distances of {:} connections'.format(inconcount, numcon))
  # print('Found {:} unused intra molecular contacts'.format(intracount))
  # print('Found {:} unused short contacts'.format(shortcount))
  # print('Found {:} unused medium molecular contacts'.format(mediumcount))
  # print('Found {:} unused long molecular contacts'.format(longcount))
  print("Finished cycle7.noa analysis")

  return assigndict2


