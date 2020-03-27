#### import the simple module from the paraview
from paraview.simple import *
import numpy as np
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

z = 0.7
R1 = 5
R2 = 4.992
x1 = np.sqrt(R1**2-z**2)
x2 = np.sqrt(R2**2-z**2)
resolution = 800
# create a new 'XML Unstructured Grid Reader'
model = 'S15_ESBC'
N = 2048
fileName = [''] * N
for i in range(0,N):
	fileName[i] = 'C:\\Users\\Zetison\\hugeFiles\\e3Dss\\paraviewResults\\'+model+'\\fluid1_time_'+str(i+1)+'.vtu'
fluid1_time_ = XMLUnstructuredGridReader(FileName=fileName)
fluid1_time_.PointArrayStatus = ['Total scalar field (real)']

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [2027, 1173]

# show data in view
fluid1_time_Display = Show(fluid1_time_, renderView1)
# trace defaults for the display properties.
fluid1_time_Display.Representation = 'Surface'
fluid1_time_Display.ColorArrayName = [None, '']
fluid1_time_Display.OSPRayScaleArray = 'Total scalar field (real)'
fluid1_time_Display.OSPRayScaleFunction = 'PiecewiseFunction'
fluid1_time_Display.SelectOrientationVectors = 'None'
fluid1_time_Display.ScaleFactor = 1.9500000000000002
fluid1_time_Display.SelectScaleArray = 'None'
fluid1_time_Display.GlyphType = 'Arrow'
fluid1_time_Display.PolarAxes = 'PolarAxesRepresentation'
fluid1_time_Display.ScalarOpacityUnitDistance = 0.3874835644455306
fluid1_time_Display.GaussianRadius = 0.9750000000000001
fluid1_time_Display.SetScaleArray = ['POINTS', 'Total scalar field (real)']
fluid1_time_Display.ScaleTransferFunction = 'PiecewiseFunction'
fluid1_time_Display.OpacityArray = ['POINTS', 'Total scalar field (real)']
fluid1_time_Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
fluid1_time_Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
fluid1_time_Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
fluid1_time_Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# create a new 'XML Unstructured Grid Reader'
fileName = [''] * N
for i in range(0,N):
	fileName[i] = 'C:\\Users\\Zetison\\hugeFiles\\e3Dss\\paraviewResults\\'+model+'\\fluid3_time_'+str(i+1)+'.vtu'
fluid3_time_ = XMLUnstructuredGridReader(FileName=fileName)
fluid3_time_.PointArrayStatus = ['Total scalar field (real)']

# show data in view
fluid3_time_Display = Show(fluid3_time_, renderView1)
# trace defaults for the display properties.
fluid3_time_Display.Representation = 'Surface'
fluid3_time_Display.ColorArrayName = [None, '']
fluid3_time_Display.OSPRayScaleArray = 'Total scalar field (real)'
fluid3_time_Display.OSPRayScaleFunction = 'PiecewiseFunction'
fluid3_time_Display.SelectOrientationVectors = 'None'
fluid3_time_Display.ScaleFactor = 0.998397643662118
fluid3_time_Display.SelectScaleArray = 'None'
fluid3_time_Display.GlyphType = 'Arrow'
fluid3_time_Display.PolarAxes = 'PolarAxesRepresentation'
fluid3_time_Display.ScalarOpacityUnitDistance = 0.32792425243476075
fluid3_time_Display.GaussianRadius = 0.499198821831059
fluid3_time_Display.SetScaleArray = ['POINTS', 'Total scalar field (real)']
fluid3_time_Display.ScaleTransferFunction = 'PiecewiseFunction'
fluid3_time_Display.OpacityArray = ['POINTS', 'Total scalar field (real)']
fluid3_time_Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
fluid3_time_Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
fluid3_time_Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
fluid3_time_Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0


# create a new 'Disk'
disk1 = Disk()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
disk1Display = Show(disk1, renderView1)
# trace defaults for the display properties.
disk1Display.Representation = 'Surface'
disk1Display.ColorArrayName = [None, '']
disk1Display.OSPRayScaleFunction = 'PiecewiseFunction'
disk1Display.SelectOrientationVectors = 'None'
disk1Display.ScaleFactor = 0.2
disk1Display.SelectScaleArray = 'None'
disk1Display.GlyphType = 'Arrow'
disk1Display.PolarAxes = 'PolarAxesRepresentation'
disk1Display.GaussianRadius = 0.1
disk1Display.SetScaleArray = [None, '']
disk1Display.ScaleTransferFunction = 'PiecewiseFunction'
disk1Display.OpacityArray = [None, '']
disk1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
disk1Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
disk1Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
disk1Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# find source
fluid3_time_ = FindSource('fluid3_time_*')

# find source
fluid1_time_ = FindSource('fluid1_time_*')

# Properties modified on disk1
disk1.OuterRadius = x1

# Properties modified on disk1
disk1.InnerRadius = x2

# Properties modified on disk1
disk1.RadialResolution = 2

# Properties modified on disk1
disk1.CircumferentialResolution = resolution

# create a new 'Calculator'
calculator1 = Calculator(Input=disk1)
calculator1.Function = ''

# Properties modified on calculator1
calculator1.Function = ''

# show data in view
calculator1Display = Show(calculator1, renderView1)
# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'
calculator1Display.ColorArrayName = [None, '']
calculator1Display.OSPRayScaleFunction = 'PiecewiseFunction'
calculator1Display.SelectOrientationVectors = 'None'
calculator1Display.SelectScaleArray = 'None'
calculator1Display.GlyphType = 'Arrow'
calculator1Display.PolarAxes = 'PolarAxesRepresentation'
calculator1Display.GaussianRadius = 0.5
calculator1Display.SetScaleArray = [None, '']
calculator1Display.ScaleTransferFunction = 'PiecewiseFunction'
calculator1Display.OpacityArray = [None, '']
calculator1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
calculator1Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
calculator1Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
calculator1Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# hide data in view
Hide(disk1, renderView1)

# Properties modified on calculator1
calculator1.CoordinateResults = 1

# Properties modified on calculator1
calculator1.Function = 'coordsX*iHat+'

# Properties modified on calculator1
calculator1.Function = 'coordsX*iHat + coordsY*jHat + (coordsZ+' + str(z) + ')*kHat'

# create a new 'Sphere'
sphere1 = Sphere()

# show data in view
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.Representation = 'Surface'
sphere1Display.ColorArrayName = [None, '']
sphere1Display.OSPRayScaleArray = 'Normals'
sphere1Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere1Display.SelectOrientationVectors = 'None'
sphere1Display.ScaleFactor = 0.1
sphere1Display.SelectScaleArray = 'None'
sphere1Display.GlyphType = 'Arrow'
sphere1Display.PolarAxes = 'PolarAxesRepresentation'
sphere1Display.GaussianRadius = 0.05
sphere1Display.SetScaleArray = [None, '']
sphere1Display.ScaleTransferFunction = 'PiecewiseFunction'
sphere1Display.OpacityArray = [None, '']
sphere1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sphere1Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sphere1Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sphere1Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# Properties modified on sphere1
sphere1.Radius = 1.0

# Properties modified on sphere1
sphere1.ThetaResolution = resolution/2

# Properties modified on sphere1
sphere1.PhiResolution = resolution

# create a new 'Sphere'
sphere2 = Sphere()

# show data in view
sphere2Display = Show(sphere2, renderView1)
# trace defaults for the display properties.
sphere2Display.Representation = 'Surface'
sphere2Display.ColorArrayName = [None, '']
sphere2Display.OSPRayScaleArray = 'Normals'
sphere2Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere2Display.SelectOrientationVectors = 'None'
sphere2Display.ScaleFactor = 0.1
sphere2Display.SelectScaleArray = 'None'
sphere2Display.GlyphType = 'Arrow'
sphere2Display.PolarAxes = 'PolarAxesRepresentation'
sphere2Display.GaussianRadius = 0.05
sphere2Display.SetScaleArray = [None, '']
sphere2Display.ScaleTransferFunction = 'PiecewiseFunction'
sphere2Display.OpacityArray = [None, '']
sphere2Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sphere2Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sphere2Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sphere2Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# Properties modified on sphere2
sphere2.Radius = R2

# Properties modified on sphere2
sphere2.ThetaResolution = resolution/2

# Properties modified on sphere2
sphere2.PhiResolution = resolution

# create a new 'Clip'
clip1 = Clip(Input=sphere2)
clip1.ClipType = 'Plane'
clip1.Scalars = [None, '']

# Properties modified on clip1
clip1.Scalars = ['POINTS', '']

# show data in view
clip1Display = Show(clip1, renderView1)
# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.OSPRayScaleArray = 'Normals'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.9899999618530274
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.6915406029087948
clip1Display.GaussianRadius = 0.4949999809265137
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# hide data in view
Hide(sphere2, renderView1)
# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 0.0, -1.0]

# Properties modified on clip1.ClipType
clip1.ClipType.Origin = [0.0, 0.0, z]

# create a new 'Sphere'
sphere3 = Sphere()

# show data in view
sphere3Display = Show(sphere3, renderView1)
# trace defaults for the display properties.
sphere3Display.Representation = 'Surface'
sphere3Display.ColorArrayName = [None, '']
sphere3Display.OSPRayScaleArray = 'Normals'
sphere3Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere3Display.SelectOrientationVectors = 'None'
sphere3Display.ScaleFactor = 0.1
sphere3Display.SelectScaleArray = 'None'
sphere3Display.GlyphType = 'Arrow'
sphere3Display.PolarAxes = 'PolarAxesRepresentation'
sphere3Display.GaussianRadius = 0.05
sphere3Display.SetScaleArray = [None, '']
sphere3Display.ScaleTransferFunction = 'PiecewiseFunction'
sphere3Display.OpacityArray = [None, '']
sphere3Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sphere3Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sphere3Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sphere3Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# Properties modified on sphere3
sphere3.Radius = R1

# Properties modified on sphere3
sphere3.ThetaResolution = resolution/2

# Properties modified on sphere3
sphere3.PhiResolution = resolution

# create a new 'Clip'
clip2 = Clip(Input=sphere3)
clip2.ClipType = 'Plane'
clip2.Scalars = [None, '']

# Properties modified on clip2
clip2.Scalars = ['POINTS', '']

# show data in view
clip2Display = Show(clip2, renderView1)
# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.OSPRayScaleArray = 'Normals'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'None'
clip2Display.SelectScaleArray = 'None'
clip2Display.GlyphType = 'Arrow'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.698525858534854
clip2Display.GaussianRadius = 0.5
clip2Display.SetScaleArray = [None, '']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = [None, '']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-0.00206091040046649, 0.0, 0.5, 0.0, -0.00019589119265828528, 0.46875, 0.5, 0.0, 0.0020624312911391206, 1.0, 0.5, 0.0]

# hide data in view
Hide(sphere3, renderView1)

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [0.0, 0.0, -1.0]

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [0.0, 0.0, z]

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# set active source
SetActiveSource(fluid1_time_)

# set scalar coloring
ColorBy(fluid1_time_Display, ('POINTS', 'Total scalar field (real)'))

# rescale color and/or opacity maps used to include current data range
fluid1_time_Display.RescaleTransferFunctionToDataRange(True, False)

# get color transfer function/color map for 'Totalscalarfieldreal'
totalscalarfieldrealLUT = GetColorTransferFunction('Totalscalarfieldreal')

# set active source
SetActiveSource(fluid3_time_)

# set scalar coloring
ColorBy(fluid3_time_Display, ('POINTS', 'Total scalar field (real)'))

# rescale color and/or opacity maps used to include current data range
fluid3_time_Display.RescaleTransferFunctionToDataRange(True, False)

# Rescale transfer function
totalscalarfieldrealLUT.RescaleTransferFunction(-1.0, 1.0)

# get opacity transfer function/opacity map for 'Totalscalarfieldreal'
totalscalarfieldrealPWF = GetOpacityTransferFunction('Totalscalarfieldreal')

# Rescale transfer function
totalscalarfieldrealPWF.RescaleTransferFunction(-1.0, 1.0)

# Properties modified on renderView1
renderView1.Background = [1.0, 1.0, 1.0]

# Properties modified on fluid3_time_Display
fluid3_time_Display.Opacity = 0.9

# set active source
SetActiveSource(fluid1_time_)

# Properties modified on fluid1_time_Display
fluid1_time_Display.Opacity = 0.9

# current camera placement for renderView1
renderView1.CameraPosition = [0.0, -18.0, 15.0]
renderView1.CameraFocalPoint = [0.0, 0.0, 0.3]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 32.8
renderView1.CameraParallelScale = 8.66025403784439
RenderAllViews()

# save screenshot
# SaveScreenshot('C:\Users\Zetison\OneDrive - SINTEF\work\movies\e3Dss\'+model+'.png', magnification=2, quality=100, view=renderView1)

# Properties modified on renderView1
renderView1.Background = [0.0, 0.0, 0.0]

# save animation images/movie
# WriteAnimation('C:\Users\Zetison\OneDrive - SINTEF\work\movies\e3Dss\'+model+'.ogv', Magnification=1, FrameRate=15.0, Compression=True)