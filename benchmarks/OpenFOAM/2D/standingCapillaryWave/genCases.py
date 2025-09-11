import casefoam


#case = [['plicRDF','mycCartesian','mlpNormal'],
#        ['RDF','heightFunction','mlpCurvature'],
#        ['mid','coarse','fine']]

case = [['mycCartesian'], #,'plicRDF','mycCartesian'],
        ['mlpCurvature','heightFunction','RDF'],
        ['coarse','coarser','coarsest','mid','fine']] #,'fine']]

#case = [['plicRDF'],
#        ['RDF','heightFunction','mlpCurvature'],
 #       ['mid','coarse','coarser','coarsest']] #,'fine']]
#        ['coarse','coarser','coarsest']] #,'fine']]


update_isoSurface = {
    'system/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'MULESScheme',
                                    'reconstructionScheme': 'isoSurface'}}}}

update_plicRDF = {
    'system/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'isoAdvection',
                                    'reconstructionScheme': 'plicRDF'}}}}

update_mycCartesian = {
    'system/fvSolution': {
        'solvers': {'alpha.water': {'advectionScheme': 'isoAdvection',
                                    'reconstructionScheme': 'mycCartesian'}}}}


update_RDF = {
    'system/curvatureDict': {
        'curvatureModel': 'RDF'}}

update_heightFunction = {
    'system/curvatureDict': {
        'curvatureModel': 'heightFunction'}}

update_interFoam = {
    'system/curvatureDict': {
        'curvatureModel': 'interFoam'}}

update_mlpCurvature = {
    'system/curvatureDict': {
        'curvatureModel': 'mlpCurvature'}}


#def updateBlockMesh(x):
#    return  {'system/blockMeshDict': {
#            'blocks': ['hex',
#                   [0, 1, 2, 3, 4, 5, 6, 7],
#                   '(%s %s 1)' % (x , 3*x),
#                   'simpleGrading',
#                   '(1 1 1)']}}

def updateBlockMesh(x,L,H,W):
    return  {'system/blockMeshDict': {
            'blocks': ['hex',
                   [0, 1, 2, 3, 4, 5, 6, 7],
                   '(%s %s %s)' % (x,1, 3*x),
                   'simpleGrading',
                       '(1 1 1)'],
        'vertices': ['(%f %f %f)' % (-L,-W,-H),
                     '(%f %f %f)' % (L,-W,-H),
                     '(%f %f %f)' % (L,W,-H),
                     '(%f %f %f)' % (-L,W,-H),
                     '(%f %f %f)' % (-L,-W,H),
                     '(%f %f %f)' % (L,-W,H),
                     '(%f %f %f)' % (L,W,H),
                     '(%f %f %f)' % (-L,W,H)]}}


data = {'plicRDF': update_plicRDF,
    'mycCartesian': update_mycCartesian,
    'RDF': update_RDF,
    'heightFunction': update_heightFunction,
    'interFoam': update_interFoam,
    'mlpCurvature': update_mlpCurvature,
    'coarsest': updateBlockMesh(8,1,3,0.25),
    'coarser': updateBlockMesh(16,1,3,0.25),
    'coarse': updateBlockMesh(32,1,3,0.25),
    'mid': updateBlockMesh(64,1,3,0.125),
    'fine': updateBlockMesh(128,1,3,0.0625)}

casefoam.mkCases('sinWave', case, data, 'tree',writeDir='Cases')

