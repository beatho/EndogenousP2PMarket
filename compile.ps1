rm build/script/cpp/*.obj
#rm build/script/cuda/*.obj
#rm build/script/cuda/lib/*
#foreach ($f in $(ls src/*.cu)){
#	$truc = $f.Name
#	nvcc -O2 --device-c  --gpu-architecture=sm_52 src/$truc -o build/script/cuda/$truc.obj "-IC:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.11_3.11.2544.0_x64__qbz5n2kfra8p0\include"
#	if ($?){
#		echo OK
#	} else{
#		echo ERROOOOOOOOOOOOOOORRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
#		exit $False
#	}
#}
nvcc -dlink $(ls build/script/cuda/*.obj)  -o build/script/cuda/lib/dlink.obj
nvcc --lib $(ls build/script/cuda/*.obj) build/script/cuda/lib/dlink.obj -o build/script/cuda/lib/dlink.lib
$str = @()

foreach ($f in $(ls src/*.cpp)){
	$truc = $f.Name
	$bidule = $truc.replace(".cpp", ".obj")
	$str += -join("build/script/cpp/",$bidule)
	cl -O2 /c /nologo /O2 /W3 /GL /DNDEBUG  -Isrc\include "-IC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.7\include" "-IC:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.11_3.11.2544.0_x64__qbz5n2kfra8p0\include" "-IC:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\include" "-IC:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\atlmfc\include" "-IC:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\VS\include" "-IC:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\ucrt" "-IC:\Program Files (x86)\Windows Kits\10\Include\10.0.19041.0\\um" "-IC:\Program Files (x86)\Windows Kits\10\\Include\10.0.19041.0\\shared" "-IC:\Program Files (x86)\Windows Kits\10\\Include\10.0.19041.0\\winrt" "-IC:\Program Files (x86)\Windows Kits\10\\Include\10.0.19041.0\\cppwinrt" /EHsc src/$truc /Fobuild/script/cpp/$($truc.replace('.cpp', '.obj')) /MT /wd4819 /wd4251 /wd4244 /wd4267 /wd4275 /wd4018 /wd4190 /EHsc
	if ($?){
		echo OK
	} else{
		echo ERROR
		exit $False
	}
	}
	$bidule=$truc.replace(".cpp", ".obj")
	$bidule=-join("build/script/cpp/", $bidule)
$str += "build/script/cuda/lib/dlink.lib"
#$str += $(ls build/script/cuda/fat/*)
echo @str
link /nologo /LTCG /INCREMENTAL:NO /DLL /MANIFEST:EMBED,ID=2 /MANIFESTUAC:NO cudart.lib cudadevrt.lib kernel32.lib user32.lib @str "/LIBPATH:C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.7\lib\x64" "/LIBPATH:C:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.11_3.11.2544.0_x64__qbz5n2kfra8p0\libs" "/LIBPATH:C:\Program Files\WindowsApps\PythonSoftwareFoundation.Python.3.11_3.11.2544.0_x64__qbz5n2kfra8p0" "/LIBPATH:C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\atlmfc\lib\x64" "/LIBPATH:C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.33.31629\lib\x64" "/LIBPATH:C:\Program Files (x86)\Windows Kits\10\lib\10.0.19041.0\ucrt\x64" "/LIBPATH:C:\Program Files (x86)\Windows Kits\10\\lib\10.0.19041.0\\um\x64" "/LIBPATH:C:\Program Files (x86)\Windows Kits\10\bin\10.0.19041.0\x64"  cudart_static.lib  /EXPORT:PyInit_EndoCuda  /OUT:build\lib.win-amd64-cpython-311\EndoCuda.cp311-win_amd64.pyd /IMPLIB:build\temp.win-amd64-cpython-311\Release\src\EndoCuda.cp311-win_amd64.lib