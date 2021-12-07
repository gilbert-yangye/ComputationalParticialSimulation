import vtk
import pytest
import numpy as np
    
def test_file_writer_output(get_num):

        reader_1 = vtk.vtkXMLPolyDataReader()
        reader_1.SetFileName("./tests/output_0.vtp")
        reader_1.Update()
        pdata_1 = reader_1.GetOutput()
        initial_num = pdata_1.GetNumberOfPoints()
        print(initial_num)
        for j in range (0,get_num):
                filename = "./tests/output_"+str(j)+".vtp"
                reader = vtk.vtkXMLPolyDataReader()
                reader.SetFileName(filename)
                reader.Update()
                pdata = reader.GetOutput()
        #assert pdata.GetNumberOfCells() == 132793
        # number of points change everytime
                assert pdata.GetNumberOfPoints() == initial_num
                for i in range(0,pdata.GetNumberOfPoints()):
                        assert pdata.GetPoint(i)[0] <=20.52  and pdata.GetPoint(i)[0] >= (-0.52) # location of points
                        assert pdata.GetPoint(i)[1] <=10.52  and pdata.GetPoint(i)[1] >= (-0.52)
                        assert pdata.GetPoint(i)[2] == 0.0
                        assert (pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0] != (-float('inf')) and \
                                pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0] != float('inf') and \
                                pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0] != np.nan)
                        assert (pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1] != (-float('inf')) and\
                                pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1] != float('inf') and\
                                pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1] != np.nan)
                        assert pdata.GetPointData().GetArray('Velocity').GetTuple(i)[2] == 0.0
                        assert (pdata.GetPointData().GetArray('Pressure').GetValue(i) != (-float('inf')) and\
                                pdata.GetPointData().GetArray('Pressure').GetValue(i) != float('inf') and\
                                pdata.GetPointData().GetArray('Pressure').GetValue(i) != np.nan)

# testing if openmp version equals to serial version
# def test_omp_output(get_num):
#         for j in range (0,get_num):
#                 filename = "output_"+str(j)+".vtp"
#                 omp_filename = "omp_output_"+str(j)+".vtp"

#                 reader = vtk.vtkXMLPolyDataReader()
#                 reader.SetFileName(filename)
#                 reader.Update()
#                 pdata = reader.GetOutput()
                
#                 omp_reader = vtk.vtkXMLPolyDataReader()
#                 omp_reader.SetFileName(omp_filename)
#                 omp_reader.Update()
#                 omp_pdata = omp_reader.GetOutput()
#     # assert pdata.GetNumberOfCells() == 132793
#     # number of points change everytime

#                 assert pdata.GetNumberOfPoints() == omp_pdata.GetNumberOfPoints()

#                 for i in range(0,pdata.GetNumberOfPoints()):
#                         assert pdata.GetPoint(i)[0] == omp_pdata.GetPoint(i)[0] 
#                         assert pdata.GetPoint(i)[1] == omp_pdata.GetPoint(i)[1]
#                         assert pdata.GetPoint(i)[2] == omp_pdata.GetPoint(i)[2]
#                         assert (pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0] == \
#                                 omp_pdata.GetPointData().GetArray('Velocity').GetTuple(i)[0])
#                         assert (pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1] == \
#                                 omp_pdata.GetPointData().GetArray('Velocity').GetTuple(i)[1])
#                         assert (pdata.GetPointData().GetArray('Velocity').GetTuple(i)[2] == \
#                                 omp_pdata.GetPointData().GetArray('Velocity').GetTuple(i)[2])
#                         assert (pdata.GetPointData().GetArray('Pressure').GetValue(i) == \
#                                 omp_pdata.GetPointData().GetArray('Pressure').GetValue(i))
                
                
@pytest.fixture()
def get_num():
        f = open('timestep.txt')             
        lines = f.readlines()                           
        num=int(lines[0])                   
        return num


