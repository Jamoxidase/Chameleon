from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SimpleLocation
from typing import Iterator

tempfile = 'tests/pSPDY.gb'

def Stealth_CDS_Extract(gbfile: str) -> Iterator:
    """
    Reads a genbank file and yields CDS/ORF features for Codon Optimization

    A generator that yields CDS features in FASTA format

    Used to produce data for codon usage statistics
    """
    gb = SeqIO.parse(gbfile,'genbank') 
    count = 0 # Used ID CDS regions, not important
    for rec in gb:
        for feature in rec.features:
            count += 1
            if feature.type == 'CDS' or feature.type == 'ORF':
                out = feature.extract(rec)

                '''Set the description to 'product' qualifier if avaiable. Else 'label' qualifier or [null] if labels are not present'''
                out.description = feature.qualifiers.get('product',feature.qualifiers.get("label",["null"]))[0]

                out.name = 'TEMP CDS HEADER'
                out.id = f'CDS_{count} CODON STATS'
                yield out
    return


def parsePlasmid(gbfile: str) -> tuple[list[SimpleLocation]]:
    avoid,cds = [],[]
    gb = SeqIO.parse(gbfile,'genbank') 
    for rec in gb:
        for feature in rec.features:
            if not (feature.type == 'CDS' or feature.type == 'ORF'):
                avoid.append((feature.location, feature.location.start % 3))
            else:
                cds.append((feature.location, feature.location.start % 3))
    return avoid,cds

one,two = parsePlasmid(tempfile)

def defineMutable(location: tuple[list[SimpleLocation]]) -> list[SimpleLocation]:
    avoid,cdsRegions = location
    temp = [] #[[[SimpLocation List],frame,parent CDS]]

    'Removes mutable cds regions from overlapping non-coding sequence'
    for cds,frame in cdsRegions:
        bounds = [[cds.start,cds.end]]
        for loc,_ in avoid:
            if cds.start > loc.start and cds.end < loc.end:
                'case: interior ORF to non-coding region'
                break
            if bounds[-1][0] < loc.start and bounds[-1][1] > loc.end:
                'case: non-coding region interior to ORF'
                new_bound = [loc.end,cds.end]
                bounds[-1][1] = loc.start
                bounds.append(new_bound)
            elif bounds[-1][0] in loc:
                'case: clipped start'
                bounds[-1][0] = loc.end
            elif bounds[-1][1] in loc:
                'case: clipped end'
                bounds[-1][1] = loc.start
        if bounds[-1][1] <= bounds[-1][0]:
            del bounds[-1]
        if not bounds:
            continue
        mut_range = [SimpleLocation(start=rng[0],end=rng[1],strand=cds.strand) for rng in bounds ]
        temp.append([mut_range,frame,cds])

    pp(temp)

    def _trimStart(loc_list: list[SimpleLocation], cds: SimpleLocation)  -> list[SimpleLocation]:
            ref = SimpleLocation(cds.start,cds.start+3,cds.strand)
            for i in range(len(loc_list)):
                loc = loc_list[i]
                if loc.start < ref.start and loc.end > ref.end:
                    loc_list[i] = SimpleLocation(loc.start,ref.start,loc.strand)
                    loc_list.insert(i+1,SimpleLocation(ref.end,loc.end,loc.strand)) 
                    break  
                if loc.start in ref:
                    loc_list[i]= SimpleLocation(ref.end,loc.end,loc.strand)
                    break
                elif loc.end in ref:
                    loc_list[i] = SimpleLocation(loc.start,ref.start,loc.strand)
                    break
            return

            
    for i in range(len(temp)):
        tLoc = temp[i][0]
        frm = temp[i][1]
        reg = temp[i][2]
        window = temp[i+1:]
        for j in range(len(window)):
            _,parent_frame,parent_cds = window[j]
            if parent_cds.start in reg:
                _trimStart(temp[i][0],parent_cds)
        trimmed = SimpleLocation(reg.start+3,reg.end,reg.strand)
        if tLoc[0].start == reg.start:
            temp[i][0][0] = SimpleLocation(trimmed.start,tLoc[0].end,tLoc[0].strand)
        temp[i][2] = trimmed
                

    'frame adjusts mutable fragments to lie in frame of parent CDS'
    def _frameAdj(loc_list: list[SimpleLocation], frame: int) -> SimpleLocation:
        'helper func'
        new = []
        for loc in loc_list:
            st = loc.start
            ed = loc.end
            for _ in range(1,3):
                if st % 3 != frame :
                    st += 1
                if ed % 3 != frame:
                    ed -= 1
            if (ed-st) < 3:
                continue
            adj = SimpleLocation(st,ed,loc.strand)
            new.append(adj)
        return new
    
    for i in range(len(temp)):
        temp[i][0] = _frameAdj(temp[i][0],temp[i][1])

    
    temp = sorted(temp,key = lambda x: len(x[2]), reverse=True)

    def _removeOverlap(loc_list: list[SimpleLocation],cds: SimpleLocation) -> list[SimpleLocation]:
        ret = []
        for loc in loc_list:
            new = [[loc.start,loc.end]]
            if loc.start in cds and loc.end in cds:
                'completely removes region'
                continue
            if loc.start < cds.start and loc.end > cds.end:
                'splits range into two'
                temp = [cds.end,loc.end]
                new[-1][1] = cds.start
                new.append(temp) 
            elif loc.start in cds:
                'truncates start'
                new[-1][0] = cds.end
            elif loc.end in cds:
                'truncates end'
                new[-1][1] = cds.start
            ret.extend([SimpleLocation(i,j,loc.strand) for i,j in new])
        return ret
    
    def _windowCorrection(old: list[list[SimpleLocation],int,SimpleLocation], st_bound: int, ed_bound: int) -> list[list[SimpleLocation],int,SimpleLocation]:
        loc,fm,parent = old
        new = SimpleLocation(st_bound,ed_bound,parent.strand)
        temp = []
        for i in loc:
            if st_bound > i.end:
                continue
            if i.start in new and i.end in new:
                temp.append(i)
            elif i.start in new:
                temp.append(SimpleLocation(i.start,ed_bound,i.strand))     
                continue
            elif i.end in new:
                temp.append(SimpleLocation(st_bound,i.end,i.strand))
                continue
            
        return [temp,fm,new]

    for i in range(len(temp)):
        frag = temp[i][0]
        frm = temp[i][1]
        region = temp[i][2]
        window = temp[i+1:]
        if region == None:
            continue
        for j in range(len(window)):
            _,parent_frame,parent_cds = window[j]
            if parent_cds == None:
                continue
            if parent_cds.start in region or parent_cds.end in region:
                if parent_frame != frm:
                   temp[i][0] = _removeOverlap(frag,parent_cds)
                if parent_frame == frm:
                    pass
            if parent_cds.start in region and parent_cds.end in region:
                temp[i+j+1][2] = None
            elif parent_cds.start in region:
                temp[i+j+1] = _windowCorrection(window[j],region.end,parent_cds.end)
            elif parent_cds.end in region:
                temp[i+j+1] = _windowCorrection(window[j],parent_cds.start,region.start)
    pp(temp)
               

            





def pp(a):
    for i,j,k in a:
        print([str(x) for x in i],k)
    print('\n')

    


    
              
           
                
defineMutable((one,two))

# print('\n\n')
# print(sorted(one,key= lambda x: len(x[0]), reverse=True))
                
                    
            