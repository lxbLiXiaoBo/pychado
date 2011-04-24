#!/usr/bin/python

#Parses an intput tabular file with: >ultracontig and then
# supercontig start, end, frame, score
#Inserts featurelocs

import sys
import psycopg2

class GoTerm(object):
    def __init__(self,id,name):
        self.id=id
        self.name=name
        self.childs=[]
        self.features=[]
        self.is_child=False

    def __repr__(self):
        return '%s (%d)' % (self.name, self.id)

    def get_all_features(self):
        f=self.features
        for c in self.childs:
            f=f+c.get_all_features()
        f=dict(map(lambda i: (i,1),f)).keys()
        return f

    def feature_count(self):
        return len(self.get_all_features())

    def add_feature(self,feature):
        self.features.append(feature)
    
    def add_child(self,child):
        self.childs.append(child)
    
    def html_dump(self,feature_link_prefix='/'):
        code="<h4>"+self.name+" (%d)</h4><ul>\n"%self.feature_count()
        for c in self.childs:
            code+="<li>"+c.html_dump(feature_link_prefix)+"</li>\n"
        for f in self.features:
            code+="<li> <a href=\"%s/%d\">%s (feature_id %d)</a></li>\n"%(feature_link_prefix, f.id, f.name,f.id)
        code+="</ul>\n"
        return code

class Feature(object):
    def __init__(self, id, name):
        self.id=id
        self.name=name

    def cv_html_from_db(self,chado,cv_prefix=''):
        chado.cur.execute("select db.name, dbx.accession, cv.name, t.name, t.definition from feature_cvterm ft, cvterm t, cv, dbxref dbx, db where feature_id=%s and ft.cvterm_id=t.cvterm_id and dbx.dbxref_id=t.dbxref_id and t.cv_id=cv.cv_id and db.db_id=dbx.db_id;" , [self.id])
        code='<ul>'
        for dbname, accession, cvname, termname, termdef in chado.cur.fetchall():
            termlink='#'
            if dbname=='GO':
                termlink='http://amigo.geneontology.org/cgi-bin/amigo/term-details.cgi?term=%s' % accession

            code+='<li> <a href="%s" target="_blank"> <h4>%s:%s</h4></a>"%s" (%s) <br/> %s </li>' % (termlink, dbname, accession, termname, cvname, termdef)
        code+='</ul>'
        return code
    
    def match_html_from_db(self,chado,cv_prefix=''):
        chado.cur.execute("select t.uniquename, t.name, md.normscore, md.rawscore from feature m, analysisfeature md, featureloc mlt, feature t, featureloc mlq where mlt.srcfeature_id=t.feature_id and mlt.feature_id=m.feature_id and mlt.rank=1 and mlq.feature_id=m.feature_id and mlq.srcfeature_id=%s and mlq.rank=0 and md.feature_id=m.feature_id order by normscore desc;" , [self.id])
        code='<table>'
        code+='<tr> <td><h4> Match accession </h4></td> <td> <h4><center>Name</center></h4> </td> <td><h4>Norm Score</h4></td> <td><h4>Raw Score</h4></td> </tr>'
        for uniquename, name, normscore, rawscore in chado.cur.fetchall():
            termlink='#'
            if uniquename.startswith('YP_') or uniquename.startswith('NP_'):
                termlink='http://www.ncbi.nlm.nih.gov/protein/%s' % uniquename
            code+='<tr> <td> <a href="%s" target="_blank"> %s</a> </td> <td>%s</td> <td>%s</td> <td>%s</td> </tr>' % (termlink, uniquename,name,normscore,rawscore)
        code+='</table>'
        return code



class PyChado(object):
    def __init__(self,database):
        self.db= psycopg2.connect(database)
        self.cur=self.db.cursor()

    def get_organisms(self):
        self.cur.execute("SELECT organism_id,common_name from organism;")
        return [(int(org[0]),org[1]) for org in self.cur.fetchall()]


    def insert_feature(self,organism,name,type_id,residues,analysis):
        self.cur.execute("insert into feature (organism_id,name,uniquename,residues,seqlen,type_id,is_analysis) values (%s,%s,%s,%s,%s,%s,%s)", [organism, name, name, residues, len(residues), type_id, analysis])

    def get_feature_id(self,name):
        self.cur.execute("SELECT feature_id from feature where name=(%s);",[name])
        return int(self.cur.fetchone()[0])
    
    def get_feature_id2(self,name):
        self.cur.execute("SELECT feature_id from feature where uniquename=(%s);",[name])
        return int(self.cur.fetchone()[0])

    def get_feature_name(self,id):
        self.cur.execute("SELECT name from feature where feature_id=%s;",[id])
        return self.cur.fetchone()[0]

    def get_feature_from_id(self, id):
        return Feature(id, self.get_feature_name(id))

    def get_cvterm_id(self,db,accession):
        self.cur.execute("select cvterm_id from db, dbxref, cvterm where db.name=%s and db.db_id=dbxref.db_id and dbxref.dbxref_id=cvterm.dbxref_id and dbxref.accession=%s;",[db,accession])
        return int(self.cur.fetchone()[0])

    def place_child(self,id, parentid, start, end, strand):
        self.cur.execute("delete from featureloc where feature_id=%s and srcfeature_id=%s;",[id, parentid])
        self.cur.execute("insert into featureloc (feature_id, srcfeature_id, fmin, fmax, strand) values (%s,%s,%s,%s,%s);",[id, parentid, start, end, strand])
        return

    def get_child_coords(self,parentid,start,end, strand):
        self.cur.execute("select feature_id,fmin,fmax,strand from featureloc where srcfeature_id=%s and fmin<=%s and fmax >=%s", [parentid,start,end])
        np=self.cur.fetchone()
        if np is None: return None
        np_id = int(np[0])
        np_start= int(np[1])
        np_end= int(np[2])
        np_strand= int(np[3])
        nc_start=start-np_start
        nc_end=end-np_start
        if strand==-1:
            nc_start=end-np_start
            nc_end=start-np_start
        if np_strand==-1:
            nc_strand=-strand
        else: nc_strand=strand
        return (np_id,nc_start,nc_end,nc_strand)

    def get_seq_supart(self,parentid,start,end,strand):
        self.cur.execute("select residues from feature where feature_id=%s", [parentid])
        res=self.cur.fetchone()[0]
        res=res[start:end]
        if strand==-1: res=res[::-1]
        return res

    def insert_feature_cvterm_file(self,file,publication=''):
        #TODO: performance
        #TODO: use ranks
        #Try to get ID for the publication, insert if not present
        self.cur.execute("select pub_id from pub where uniquename=%s;",[publication])
        try:
            pub_id=self.cur.fetchone()[0]
        except:
            self.cur.execute("insert into pub (uniquename,type_id) values (%s,1);", [publication])
            self.db.commit();
            self.cur.execute("select pub_id from pub where uniquename=%s;",[publication])
            pub_id=self.cur.fetchone()[0]
        for line in file.readlines():
            try:
                if line[-1]=='\n': line=line[:-1]
            except:
                pass
            line=line.split(',')
            dbname,accession=line[1].split(':')
            #get ids
            feature_id= get_feature_id(line[0])
            cvterm_id= get_cvterm_id(dbname,accession)
            self.cur.execute("insert into feature_cvterm (feature_id,cvterm_id,pub_id) values (%s,%s,%s);", [feature_id,cvterm_id,pub_id])
        self.db.commit()

    def insert_blast_xml_results(self,organism_id,file,analysis='blast run',description=''):
        from Bio.Blast import NCBIXML
        somatch_id=self.get_id_cv_term('sequence','match')
        somatchpart_id=self.get_id_cv_term('sequence','match_part')
        region_id=self.get_id_cv_term('sequence','region')

        analysisid=None
        for query in NCBIXML.parse(file):
            if analysisid is None:
                #Create a record for the blast run (analysis)
                self.cur.execute('insert into analysis (name,description,program,programversion) values (%s,%s,%s,%s)', [analysis,description,query.application,query.version])
                self.cur.execute('select analysis_id from analysis where name=%s order by analysis_id desc limit 1',[analysis]);
                analysisid=self.cur.fetchone()[0]

            if len(query.descriptions)>0:
                #check for the query to be on the database
                qfeatureid=self.get_feature_id(query.query.split()[0])
                assert qfeatureid > 0
            
            for i in xrange(len(query.descriptions)):
                #insert feature for the hit from description SO:match
                featurename='blast_match: %s to %s from %s (%d)'%(query.query.split()[0],query.descriptions[i].title[:100],analysis,analysisid)
                featurename=featurename[:255]
                self.cur.execute('insert into feature (organism_id,uniquename,name,is_analysis,type_id) values (%s,%s,%s,%s,%s)',
                        [organism_id,featurename,query.descriptions[i].title[:100],True,somatch_id]);
                featureid=self.get_feature_id2(featurename)
                
                #insert featureanalysis with scores from desc.bits desc.score TODO:(where to get e-value and frac_identical)?
                self.cur.execute('insert into analysisfeature (feature_id,analysis_id,rawscore,normscore) values (%s,%s,%s,%s)',
                    [featureid,analysisid,query.descriptions[i].bits,query.descriptions[i].score])
                self.cur.execute('select analysisfeature_id from analysisfeature where feature_id=%s and analysis_id=%s',[featureid,analysisid])
                analysisfeatureid=self.cur.fetchone()[0];

                #insert featureloc to query (rank=0)
                self.cur.execute('insert into featureloc (feature_id,srcfeature_id) values (%s,%s)', [featureid,qfeatureid])
                #check and insert target feature
                self.cur.execute('select feature_id from feature where uniquename=%s',[query.alignments[i].accession])
                try:
                    targetid=self.cur.fetchone()[0];
                except:
                    self.cur.execute('insert into feature (organism_id,uniquename,name,type_id) values (%s,%s,%s,%s)',[organism_id,query.alignments[i].accession,query.alignments[i].hit_def,region_id])
                    self.cur.execute('select feature_id from feature where uniquename=%s',[query.alignments[i].accession])
                    targetid=self.cur.fetchone()[0];


                #insert featureloc to target (rank=1)
                self.cur.execute('insert into featureloc (feature_id,srcfeature_id,rank) values (%s,%s,1)', [featureid,targetid])

                #for hsp in query.hsp:
                    #create the hsp feature SO:match_part
                    #create the hsp featureanalysis
                    #insert featureloc to query (rank=0)
                    #insert featureloc to target (rank=1)
        self.db.commit()



        
    def get_id_cv_term(self, cv, term):
        self.cur.execute("select cvterm_id from cv v,cvterm t where v.cv_id=t.cv_id and v.name=%s and t.name=%s;",[cv,term])
        return self.cur.fetchone()[0]

    def create_cv_tree(self,cvname):
        self.cur.execute("select cv_id from cv where name=%s;",[cvname])
        cv_id=self.cur.fetchone()[0]
        
        self.cur.execute("select cvterm_id from cv v,cvterm t where v.cv_id=t.cv_id and v.name='relationship' and t.name='is_a';")
        isa_id=self.cur.fetchone()[0]
        
        terms={}
        features={}

        #get all features annotated with the cv, insert them into features and the direct annotation in terms
        self.cur.execute("select f.feature_id, f.name, t.cvterm_id, t.name from feature f, feature_cvterm ft, cvterm t where f.feature_id=ft.feature_id and t.cvterm_id=ft.cvterm_id and t.cv_id=%s" , [cv_id])
        for fid, fname, tid, tname in self.cur.fetchall():
            if not features.has_key(fid):
                features[fid]=Feature(fid,fname)
            if not terms.has_key(tid):
                terms[tid]=GoTerm(tid, tname)

            terms[tid].add_feature(features[fid])

        #recreate term tree (XXX: no cycle detection!!!)
        self.cur.execute("select r.subject_id, st.name, r.object_id, ot.name from cvterm_relationship r, cvterm st, cvterm ot where r.type_id=%s and st.cvterm_id=r.subject_id and ot.cvterm_id=r.object_id and st.cv_id=%s and ot.cv_id=%s;" , [isa_id, cv_id, cv_id])
        rels= self.cur.fetchall()
        parents_added=True
        while parents_added:
            parents_added=False
            newchilds=[]
            for sid, sname, oid, oname in rels:
                #only insert if child exists and hasn't got a father
                if terms.has_key(sid) and terms[sid].is_child==False:
                    #create father if needed
                    if not terms.has_key(oid):
                        parents_added=True
                        terms[oid]=GoTerm(oid, oname)
                    terms[oid].add_child(terms[sid])
                    newchilds.append(terms[sid])
            for c in newchilds: c.is_child=True

        return [x for x in terms.values() if not x.is_child]








