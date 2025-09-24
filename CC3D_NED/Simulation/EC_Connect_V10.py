
from cc3d import CompuCellSetup

# from NDJ_CellMorph import NDJ_CellMorph
from NDJ_SBML_Step import NDJ_SBML_Step
from NDJ_Interactions import NDJ_Interactions
from NDJ_Vis_DataTrack import NDJ_Vis_DataTrack

# CompuCellSetup.register_steppable(steppable=NDJ_CellMorph(frequency=1))
CompuCellSetup.register_steppable(steppable=NDJ_SBML_Step(frequency=1))
CompuCellSetup.register_steppable(steppable=NDJ_Interactions(frequency=1))
CompuCellSetup.register_steppable(steppable=NDJ_Vis_DataTrack(frequency=1))

CompuCellSetup.run()
