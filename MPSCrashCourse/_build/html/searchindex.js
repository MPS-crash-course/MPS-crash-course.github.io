Search.setIndex({"alltitles": {"Bringing It All Together": [[7, null]], "Canonical Form and Truncation": [[3, null]], "Code: Basic Exact Diagonalization": [[2, null]], "Code: Extend the MPS Class": [[1, null], [1, null]], "Code: MPS Class": [[1, null]], "Code: Product State": [[1, null]], "Density Matrix Renormalization Group (DMRG)": [[0, "density-matrix-renormalization-group-dmrg"], [5, null]], "Exercise: Contracting tensors": [[1, null]], "Experimental motivation": [[2, "experimental-motivation"]], "Introduction": [[0, null]], "MPS Crash Course": [[8, null]], "MPS from a state vector": [[1, "mps-from-a-state-vector"]], "MPS to vector": [[1, "mps-to-vector"]], "Matrix Product Operators": [[6, null]], "Matrix Product States": [[1, null], [1, "id1"]], "Matrix Product States (MPS)": [[0, "matrix-product-states-mps"]], "Numpy tensordot and transpose": [[1, "numpy-tensordot-and-transpose"]], "Product states": [[1, "product-states"]], "Quantum states of many spins": [[1, "quantum-states-of-many-spins"]], "References": [[0, "references"], [2, "references"]], "Summary": [[1, "summary"]], "Tensor Networks (TNs)": [[0, "tensor-networks-tns"]], "Tensor network diagrams": [[1, "tensor-network-diagrams"]], "Tests: MPS basics": [[1, null]], "The Heisenberg Antiferromagnet": [[2, null]], "The Heisenberg model": [[2, "the-heisenberg-model"]], "The dynamical structure factor": [[2, "the-dynamical-structure-factor"]], "Time Evolving Block Decimation (TEBD)": [[0, "time-evolving-block-decimation-tebd"], [4, null]], "Week 1": [[8, null]], "Week 2": [[8, null]], "Week 3": [[8, null]], "Week 4": [[8, null]]}, "docnames": ["content/week1/intro", "content/week1/mps", "content/week1/problem", "content/week2/canonical", "content/week2/tebd", "content/week3/dmrg", "content/week3/mpo", "content/week4/dsf", "home"], "envversion": {"sphinx": 62, "sphinx.domains.c": 3, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 9, "sphinx.domains.index": 1, "sphinx.domains.javascript": 3, "sphinx.domains.math": 2, "sphinx.domains.python": 4, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx.ext.intersphinx": 1, "sphinxcontrib.bibtex": 9}, "filenames": ["content/week1/intro.md", "content/week1/mps.md", "content/week1/problem.md", "content/week2/canonical.md", "content/week2/tebd.md", "content/week3/dmrg.md", "content/week3/mpo.md", "content/week4/dsf.md", "home.md"], "indexentries": {}, "objects": {}, "objnames": {}, "objtypes": {}, "terms": {"": [0, 1, 2], "0": [1, 2], "00": 1, "000": 1, "0000": 1, "001": 1, "01": 1, "0101": 1, "012": 0, "013": 0, "040502": 0, "045003": 0, "06": 0, "09": 0, "0th": 1, "1": [0, 1, 2], "10": [0, 1, 2], "100": [0, 1], "1007": 0, "1016": 0, "11": 1, "1103": [0, 2], "111": 2, "115": 0, "117": 0, "13": 2, "137205": 2, "147902": 0, "158": 0, "192": 0, "1988": 0, "1992": 0, "1d": [1, 2], "1j": 2, "1st": 1, "2": [0, 1, 2], "2003": 0, "2004": 0, "2010": 0, "2011": 0, "2013": 2, "2014": 0, "2021": 0, "281": 0, "2863": 0, "2866": 0, "287": 0, "3": [0, 1], "326": 0, "349": 0, "4": [1, 2], "477": 0, "5": 1, "528": 0, "5th": 1, "6": 1, "69": 0, "7": 1, "8": 1, "87": 0, "91": 0, "93": 0, "96": 0, "A": [0, 1, 2], "As": 1, "At": 8, "By": [0, 2], "For": [1, 2], "If": [0, 1, 8], "In": [0, 1, 2, 8], "It": [1, 8], "One": 1, "That": 1, "The": [0, 1, 8], "There": [2, 8], "These": [1, 2], "To": [1, 2], "_": [1, 2], "_3": 2, "__init__": 1, "_i": 2, "a_": 1, "abl": [0, 1, 2], "about": [0, 1, 2], "abov": [1, 2], "abstract": [0, 1], "access": 0, "accur": 2, "achiev": 2, "act": 1, "actual": [1, 2], "ad": [0, 1], "add": 1, "affleck": 0, "ag": 0, "again": 1, "against": 2, "agreement": 2, "aklt": 0, "aklt88": 0, "al": 2, "algorithm": [0, 1, 2, 8], "align": 2, "all": [0, 1, 2, 8], "allclos": 1, "allow": [0, 1, 2, 8], "along": [1, 8], "alongsid": 0, "alpha": 2, "alpha_": 1, "alpha_0": 1, "alpha_1": 1, "alpha_2": 1, "alpha_3": 1, "alpha_4": 1, "alpha_i": 2, "alpha_j": 2, "alpha_n": 1, "alreadi": 1, "also": [0, 1, 2], "although": 8, "am": 8, "amplitud": [0, 1], "an": [0, 1, 2], "ani": 8, "annal": 0, "ansatz": [0, 2], "anti": 2, "antiferromagnet": [0, 8], "anyth": 1, "aop": 0, "ap": [0, 2], "appear": 1, "append": 1, "appli": [0, 2], "applic": 1, "approach": [0, 1], "approxim": [0, 1], "ar": [0, 1, 2, 8], "area": 0, "argument": 1, "around": [1, 8], "arrai": [1, 2], "arriv": 1, "articl": 0, "aspect": [2, 8], "assert": 1, "associ": 1, "assum": 8, "attempt": 1, "attribut": 1, "audienc": 8, "ax": 1, "b": [0, 1, 2], "b_": 1, "back": 1, "barthel": 2, "base": [0, 8], "basi": 1, "basic": 8, "becaus": [1, 2], "becom": [0, 1], "befor": [0, 1, 2], "behind": 0, "below": 1, "beta": 2, "beta_k": 2, "beth": 2, "better": 1, "between": [0, 1], "beyond": 2, "bf01218021": 0, "bf01309281": 0, "binari": 1, "bit": 0, "black": 8, "block": [1, 2, 8], "bodi": [0, 1, 8], "bond": [0, 1], "both": [1, 2], "box": 8, "branch": 8, "bring": 8, "broken": 2, "built": 0, "c": [1, 2], "c_": 1, "calcul": [1, 2], "call": [1, 2, 8], "came": 0, "can": [0, 1, 2, 8], "candid": 2, "canon": [1, 8], "case": [0, 1, 2], "caux": 2, "cdot": [1, 2], "centr": [0, 1], "central": 1, "certain": 2, "chain": [0, 1, 2], "channel": 0, "check": [1, 2], "chemistri": 0, "chi": 1, "chi_": 1, "chi_n": 1, "choic": 1, "choos": 1, "cirac": 0, "circl": 1, "circuit": 0, "cl": 1, "classic": [0, 1], "classmethod": 1, "clearer": 0, "code": 8, "coeffici": 1, "collect": 1, "colour": 1, "column": 1, "com": 0, "combin": [1, 2], "come": 1, "comment": 1, "commonli": 0, "commun": 0, "compar": 2, "complet": 1, "complex": 1, "compon": 2, "compound": 2, "comput": [0, 1, 2], "concept": 0, "concern": 2, "conda": 8, "condens": 0, "condit": 1, "confus": 1, "conjug": 1, "connect": 1, "consid": [0, 1, 2], "consist": 1, "constant": 2, "constantli": 1, "construct": [1, 2], "content": 8, "context": [0, 1], "continu": 1, "continua": 2, "continuum": 2, "contract": 0, "control": 0, "convert": 1, "copi": 1, "core": 1, "correctli": 1, "correl": 2, "correspond": 1, "cost": 0, "costli": 2, "could": 1, "coupl": 2, "cours": [0, 1, 2], "cover": [2, 8], "cperezgarciasv21": 0, "creat": [1, 8], "cross": 2, "crucial": 2, "cumbersom": 1, "curs": 1, "d": 2, "dagger": 1, "dangl": 1, "dash": 1, "data": 2, "david": 0, "deal": [0, 1], "dealt": 1, "dec": 0, "decim": [1, 2, 8], "decompos": 1, "decomposit": 1, "def": [1, 2], "defin": [1, 2], "definit": 1, "deflect": 2, "degre": 1, "demystifi": [0, 8], "densiti": [2, 8], "describ": 2, "desir": 0, "detail": [0, 1, 2], "develop": 8, "diagon": [0, 1], "diagramat": 1, "did": 0, "differ": [0, 1], "dimens": [0, 1, 2], "dimension": [0, 1, 2], "directori": 1, "discuss": [0, 1, 2], "dispers": 2, "dispos": 1, "distinct": 2, "distribut": 0, "dmrg": [1, 2, 8], "do": [0, 1, 2], "doi": [0, 2], "done": [1, 2], "down": 1, "downarrow": 1, "draw": 1, "drawn": 1, "dssf": 2, "dt": 2, "due": [0, 2], "dure": 0, "dynam": 0, "e": [1, 2], "each": [0, 1, 8], "easiest": 1, "easili": [1, 2], "ed": [1, 2], "edg": 0, "effect": 1, "effici": [0, 1], "eigsh": 2, "electron": 1, "element": 1, "elliott": 0, "els": 1, "encourag": 1, "end": [1, 2], "energi": [0, 2], "entangl": 0, "entropi": 0, "environ": 8, "eq": 1, "equal": 1, "equat": 1, "equival": 1, "error": 0, "et": 2, "ever": 0, "everyon": 8, "everyth": 8, "evolut": [0, 1, 2, 8], "evolv": [1, 2, 8], "exact": [0, 1], "exactli": [0, 1, 2], "exampl": 1, "except": 1, "excit": 2, "exclud": 1, "expect": [0, 2], "experi": 2, "expert": 8, "explicitli": 1, "exponenti": [1, 2], "expos": 1, "express": [0, 1], "extend": 2, "extrem": [0, 2], "ey": 2, "fact": 2, "factor": 1, "fail": 1, "fals": 1, "familiar": 1, "far": 0, "fashion": 1, "feedback": 8, "feel": [1, 8], "fig": [0, 1, 2], "figur": 1, "file": [1, 2, 8], "final": 1, "find": [0, 1, 8], "finit": [0, 1, 2], "first": [0, 1, 8], "fix_path": 1, "focu": 0, "focuss": 1, "folder": 8, "follow": [0, 1, 8], "form": [0, 1, 2, 8], "formul": 0, "found": 1, "four": 8, "fourier": 2, "frac": 2, "frank": 0, "free": 1, "freedom": 1, "from": [0, 2], "fromvector": 1, "frost": 2, "full": 0, "full_matric": 1, "function": [1, 2], "fundament": 1, "furthermor": 8, "futur": 8, "f\u00fcr": 0, "g": 1, "gain": 1, "gapless": 2, "garc": 0, "gate": 0, "gener": [0, 1, 2, 8], "generalis": 1, "get": [0, 1, 8], "github": [1, 8], "given": [1, 2], "go": [0, 2, 8], "goal": [2, 8], "goe": 1, "good": [0, 1, 2], "graph": 0, "graphic": 1, "grate": 8, "ground": [0, 1, 2, 8], "groundstat": 0, "group": [1, 2, 8], "grow": 1, "growth": 2, "guid": 2, "guifr\u00e9": 0, "h": [0, 2], "ha": [0, 1, 2], "had": 1, "hal": 0, "half": 0, "halv": 0, "hamiltonian": 2, "hand": 1, "have": [1, 2], "haven": 1, "heisenberg": [1, 8], "heisenberggroundst": 2, "heisenberghamiltonian": 2, "helpfulli": 8, "henc": 1, "here": 1, "high": [2, 8], "higher": [0, 1], "highest": 0, "hilbert": [0, 2], "histor": 0, "hope": 8, "hopefulli": 0, "horizont": 1, "how": [0, 1, 8], "howev": [0, 1, 2, 8], "http": [0, 2], "i": [0, 1, 2, 8], "i_1": 1, "i_2": 1, "i_3": 1, "i_n": 1, "ian": 0, "idea": [0, 1], "ideal": 2, "ignacio": 0, "iht": 2, "ij": 1, "ijk": 1, "ijkl": 1, "ik": 1, "imaginari": 0, "immedi": 1, "implement": 1, "import": [0, 1, 2], "importantli": 0, "improv": 8, "includ": [1, 2, 8], "increas": 0, "inde": 0, "index": 1, "indic": 1, "induc": 0, "inelast": 2, "infeas": 1, "infinit": 0, "inform": 2, "infti": 2, "ingredi": 0, "insight": 0, "instanc": 1, "instead": [0, 1], "int": 1, "int_": 2, "intend": 8, "interact": [1, 2], "interchang": 1, "interest": 0, "intract": 0, "intrins": 0, "introduc": [0, 1, 2], "introduct": 8, "involv": 1, "iq": 2, "isotrop": [0, 2], "issu": 0, "itensor": 8, "iter": 8, "itself": 1, "j": [0, 1, 2], "januari": 0, "jargon": 0, "ji": 1, "jul": 0, "julia": 8, "just": [1, 2], "k": [1, 2], "kcuf": 2, "keep": [1, 2], "kei": [0, 1], "kennedi": 0, "kept": 0, "kind": 1, "kj": 1, "klumpersz92": 0, "kl\u00fcmper": 0, "know": 1, "known": 1, "kron": 2, "l": [1, 2], "l1": 1, "l2": 1, "la": 1, "label": 1, "lake": 2, "lambda_i": 1, "landmark": 2, "langl": [1, 2], "languag": 0, "laptop": 1, "larg": [0, 1], "largest": [0, 1], "ldot": 1, "lead": 1, "learn": 2, "lectur": 8, "led": 0, "left": [1, 2], "leg": 1, "length": [1, 2], "let": [0, 1, 2], "lett": [0, 2], "level": [0, 1, 2], "lieb": 0, "like": [1, 8], "linalg": [1, 2], "line": 1, "link": [0, 2], "list": 1, "local": [0, 2], "log2": 1, "look": 0, "low": [1, 2], "ltc": 2, "m": 1, "m_": 1, "magnet": 2, "magnon": 2, "mai": 1, "main": [0, 1, 2, 8], "mainli": 0, "make": [0, 1], "mani": [0, 8], "manipul": 1, "match": 1, "materi": 2, "mathbf": 2, "mathemat": 0, "matric": [0, 1, 2], "matrix": [2, 8], "matter": 0, "maximis": 0, "me": 0, "mean": 1, "measur": 2, "mention": 1, "mera": 0, "method": [0, 1, 2], "middl": 1, "min": 1, "mind": 2, "minim": [1, 8], "minimum": 8, "mod": 0, "model": [0, 1, 8], "modern": 0, "momentum": 2, "more": [0, 1, 2, 8], "most": [0, 1, 8], "mp": 2, "mps_basic": 1, "mps_diagram": [], "mult": [], "multi": [0, 2], "multipl": 1, "multipli": 1, "multispinon": 2, "must": [1, 8], "n": [1, 2], "nagler": 2, "name": [0, 2], "natur": 0, "ndim": 1, "nearest": 2, "necessari": 2, "need": 1, "neg": 1, "neighbour": 2, "network": 8, "neutron": 2, "new": 1, "next": 1, "non": 1, "norbert": 0, "normal": 1, "normalis": 1, "notat": 1, "note": 1, "nov": 0, "now": [0, 1, 2], "np": [1, 2], "number": [0, 1], "numer": [0, 1, 2], "numpi": 2, "object": [0, 1], "obtain": 2, "obviou": 1, "oct": 0, "off": [1, 2], "often": 1, "omega": 2, "one": [0, 1, 2], "onli": [0, 1, 2], "oper": [0, 2, 8], "opposit": 1, "order": 1, "org": [0, 2], "origin": [0, 1], "oru14": 0, "or\u00fa": 0, "other": [0, 1], "otim": 1, "our": [1, 2], "out": [1, 8], "outcom": 1, "outgo": 2, "over": [1, 8], "own": [2, 8], "p": 1, "packag": 8, "pair": 0, "panel": 2, "paper": 0, "paradigmat": 2, "paramet": 1, "part": [1, 2], "particl": [1, 2], "particular": [0, 1, 2], "particularli": 1, "pauli": 2, "pep": 0, "perform": [0, 1, 2, 8], "perhap": 1, "permut": 1, "perturb": 2, "phy": [0, 2], "physic": [0, 1], "physik": 0, "physrevlett": [0, 2], "pi": 2, "pii": 0, "plai": 1, "point": [1, 2], "polar": 2, "posit": 1, "possibl": 1, "potenti": 0, "power": [0, 1], "practic": [0, 1, 8], "predict": 2, "prefer": 2, "present": 1, "previou": 1, "print": 1, "probabl": [0, 1], "problem": 1, "process": [0, 1], "produc": [2, 8], "product": [2, 8], "productst": 1, "prohibit": 2, "project": 0, "proof": 1, "properti": 0, "provabl": 0, "provid": [1, 2, 8], "psi": 1, "psi_": 1, "psi_0": [1, 2], "psi_1": 1, "psi_2": 1, "psi_i": 1, "psi_n": 1, "pure": 1, "put": 1, "py": [1, 2], "python": [1, 8], "p\u00e9rez": 0, "q": 2, "quad": 1, "quantiti": 2, "quantum": [0, 2, 8], "quench": 1, "quick": 1, "quickli": 1, "quit": 1, "r": [0, 1, 2], "r1": 1, "r2": 1, "randn": 1, "random": 1, "rang": 2, "rangl": [1, 2], "rank": 1, "rather": 1, "real": [0, 2], "realiz": 0, "realli": 8, "reason": 1, "reault": 1, "recommend": [0, 1, 8], "rectangular": 1, "reduc": [0, 1], "refer": 1, "regular": 1, "relat": [0, 1, 2], "remain": 1, "remark": 2, "renorm": [2, 8], "repeat": [0, 1], "repeatedli": 1, "repositori": [1, 8], "repres": [0, 1, 8], "represent": [0, 1], "requir": 1, "research": 8, "reshap": 1, "rest": 1, "restrict": [0, 1], "result": [0, 1, 2], "result1": 1, "result2": 1, "return": [1, 2], "rev": [0, 2], "reveal": 2, "review": 0, "revmodphi": 0, "rewritten": 1, "right": [1, 2], "rom\u00e1n": 0, "root_dir": 1, "row": 1, "run": 8, "s0003491610001752": 0, "s0003491614001596": 0, "s_": 2, "s_i": 2, "s_x": 2, "s_y": 2, "s_z": 2, "sa": 2, "said": 2, "same": [0, 1], "sampl": 2, "satisfi": 1, "scalar": 1, "scale": 0, "scatter": 2, "sch11": 0, "schadschneid": 0, "schemat": 0, "schmidt": 1, "schollw\u00f6ck": [0, 2], "schuch": 0, "scienc": 0, "sciencedirect": 0, "scipi": [1, 2], "scope": 2, "second": [0, 1], "section": [1, 2], "self": 1, "sens": 0, "sep": 2, "separ": 1, "seri": 0, "set": 8, "sever": 0, "shape": 1, "shorthand": 1, "should": [0, 1, 2], "show": [1, 2], "shown": [0, 1, 2], "sigma": 2, "similar": [1, 2], "simpl": [1, 2], "simpler": 1, "simpli": 1, "simplic": 2, "simplifi": [1, 2], "simul": [0, 1], "sinc": [0, 1, 2], "singl": [1, 2], "singular": 1, "site": [0, 1, 2], "size": [0, 1], "slightli": 0, "smallest": 1, "so": [1, 8], "solut": 1, "solv": [1, 2], "solvabl": 2, "some": [0, 1, 2], "soon": 1, "sp": 2, "space": [0, 2, 3, 4, 5, 6, 7], "spars": 2, "speak": 1, "speaker": 8, "special": 0, "specif": [0, 1, 2], "specifi": 1, "spectrum": 0, "spin": 2, "split": [1, 2], "squar": 1, "src": [1, 2], "start": [0, 1, 2, 8], "state": [2, 8], "step": 1, "steven": 0, "store": 1, "stress": 8, "strictli": 1, "strongli": 8, "structur": [0, 1, 8], "studi": [0, 1, 2], "subspac": 0, "success": [0, 1], "sum": 1, "sum_": [1, 2], "sum_i": 1, "sum_k": 1, "sum_p": 1, "superposit": 1, "superscript": 1, "sure": 1, "svd": 1, "symmetri": [0, 2], "system": [0, 1, 2], "t": [0, 1, 2, 8], "take": 1, "taken": 0, "talk": 8, "tasaki": 0, "task": 2, "team": 2, "tebd": [1, 2, 8], "technic": 1, "techniqu": 2, "temperatur": 2, "templat": 8, "tennant": 2, "tenpi": 8, "tensor": 8, "term": [0, 1], "text": [1, 2], "th": 1, "than": 1, "thei": [0, 1], "them": 1, "theorem": 0, "theoret": 2, "theori": 2, "therefor": [0, 1], "thermal": 0, "theta": 1, "theta_": 1, "thi": [0, 1, 2, 3, 4, 5, 6, 7, 8], "thing": 0, "think": 1, "third": 1, "those": [1, 8], "three": 1, "through": 1, "throughout": [1, 2], "time": [1, 2, 8], "togeth": 8, "tom": 0, "tool": 1, "topic": 0, "touch": 2, "tovector": 1, "toward": 1, "tr": 1, "track": 1, "tractabl": 0, "transform": 2, "treat": 1, "tree": 0, "truncat": [0, 1, 8], "try": [1, 8], "ttn": 0, "tupl": 1, "two": [0, 1, 2], "type": [0, 1], "typic": 1, "u": [0, 1, 2], "ulrich": 0, "under": 2, "underli": 0, "understand": 1, "unequ": 2, "uniform": 2, "unitari": [1, 2, 8], "unless": 8, "until": 1, "up": [0, 8], "uparrow": 1, "url": [0, 2], "us": [0, 1, 2, 8], "v": [1, 2], "v_i": 1, "valenc": 0, "valu": [0, 1], "variat": 0, "vb": 0, "vdg": 1, "veri": [0, 1, 2, 8], "version": 1, "verstraet": 0, "vertex": 0, "vertic": 1, "vid03": 0, "vid04": 0, "vidal": 0, "view": 0, "virtual": 1, "wa": 0, "wai": [0, 1], "want": [1, 2], "watch": [3, 4, 5, 6, 7], "wave": 2, "wavefunct": 2, "we": [0, 1, 2, 8], "week": [0, 1, 2], "week0": 8, "well": [0, 1, 8], "were": 0, "what": [0, 1], "when": [0, 1], "where": [0, 1, 2, 8], "wherea": 2, "whi92": 0, "which": [0, 1, 2], "while": [0, 2], "white": 0, "widetild": 1, "without": 8, "won": [0, 8], "work": [1, 2], "would": [0, 1, 8], "write": [1, 8], "written": 1, "www": 0, "x": 2, "y": 2, "yaml": 8, "you": [0, 1, 2, 8], "your": [1, 2, 8], "yourself": 1, "z": 2, "z_": 2, "zeitschrift": 0, "zero": [1, 2], "zittartz": 0, "zz": 2, "\u0131": 0}, "titles": ["<span class=\"section-number\">1. </span>Introduction", "<span class=\"section-number\">3. </span>Matrix Product States", "<span class=\"section-number\">2. </span>The Heisenberg Antiferromagnet", "<span class=\"section-number\">4. </span>Canonical Form and Truncation", "<span class=\"section-number\">5. </span>Time Evolving Block Decimation (TEBD)", "<span class=\"section-number\">7. </span>Density Matrix Renormalization Group (DMRG)", "<span class=\"section-number\">6. </span>Matrix Product Operators", "<span class=\"section-number\">8. </span>Bringing It All Together", "MPS Crash Course"], "titleterms": {"1": 8, "2": 8, "3": 8, "4": 8, "It": 7, "The": 2, "all": 7, "antiferromagnet": 2, "basic": [1, 2], "block": [0, 4], "bring": 7, "canon": 3, "class": 1, "code": [1, 2], "contract": 1, "cours": 8, "crash": 8, "decim": [0, 4], "densiti": [0, 5], "diagon": 2, "diagram": 1, "dmrg": [0, 5], "dynam": 2, "evolv": [0, 4], "exact": 2, "exercis": 1, "experiment": 2, "extend": 1, "factor": 2, "form": 3, "from": 1, "group": [0, 5], "heisenberg": 2, "introduct": 0, "mani": 1, "matrix": [0, 1, 5, 6], "model": 2, "motiv": 2, "mp": [0, 1, 8], "network": [0, 1], "numpi": 1, "oper": 6, "product": [0, 1, 6], "quantum": 1, "refer": [0, 2], "renorm": [0, 5], "spin": 1, "state": [0, 1], "structur": 2, "summari": 1, "tebd": [0, 4], "tensor": [0, 1], "tensordot": 1, "test": 1, "time": [0, 4], "tn": 0, "togeth": 7, "transpos": 1, "truncat": 3, "vector": 1, "week": 8}})