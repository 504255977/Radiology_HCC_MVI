import shutil
import subprocess
import sys
from os import listdir
from os.path import join
import os

def is_image_file(filename):
    return any(filename.endswith(extension) for extension in [".nii.gz"])

def AffBsp(mov_aff, out_aff, move_mask):
    if not os.path.exists(out_aff):
        os.makedirs(out_aff)
    cmd = sys.path[0]+"\elastix\elastix.exe -f " + fix + ' -fMask ' + fixmask + " -m " + mov_aff + ' -mMask ' + move_mask + " -out " + out_aff + " -p "+sys.path[0]+"\elastix\\Par0057Bspline.txt"
    subprocess.call(cmd, shell=True)

    # if not os.path.exists(out_bsp):
    #     os.makedirs(out_bsp)
    # cmd = sys.path[0]+"\elastix\elastix.exe -f " + fix + " -m " + mov_bsp + " -out " + out_bsp + " -p "+sys.path[0]+"\elastix\\Par0047Pairwise.txt"
    # subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    # 读取图像的路径
    input_time0_dir = 'input_data/aAif'
    input_time1_dir = 'input_data/dAif'
    input_time2_dir = 'input_data/nAif'
    input_time3_dir = 'input_data/v'
    # 读取mask的路径
    input_mask0_dir = 'input_data_mask/aAif_mask'
    input_mask1_dir = 'input_data_mask/dAif_mask'
    input_mask2_dir = 'input_data_mask/nAif_mask'
    input_mask3_dir = 'input_data_mask/v_mask'
    # 配准结果的存放路径
    output_time0_dir = 'output_data/aRes'
    output_time1_dir = 'output_data/dRes'
    output_time2_dir = 'output_data/nRes'

    tmp_time0_dir = 'tmp/a'
    tmp_time1_dir = 'tmp/d'
    tmp_time2_dir = 'tmp/n'

    time3_filenames = [x for x in listdir(input_time3_dir)if is_image_file(x)]
    mask3_filenames = [x for x in listdir(input_mask3_dir)if is_image_file(x)]
    num = len(time3_filenames)
    for i in range(num):
        fix = join(input_time3_dir, time3_filenames[i])
        fixmask = join(input_mask3_dir, mask3_filenames[i])

        # mask-->time1
        mov_aff = join(input_time0_dir, time3_filenames[i])
        out_aff = join(tmp_time0_dir, time3_filenames[i][:-7]+'_aff')
        move_mask = join(input_mask0_dir, mask3_filenames[i])
        # mov_bsp = join(out_aff, 'result.0.nrrd')
        # out_bsp = join(tmp_mask_dir, filename[:-5]+'_bsp')
        AffBsp(mov_aff, out_aff, move_mask)

        shutil.move(join(out_aff, 'result.0.nii'), join(output_time0_dir, time3_filenames[i]))

        # time2-->time1
        mov_aff = join(input_time1_dir, time3_filenames[i])
        out_aff = join(tmp_time1_dir, time3_filenames[i][:-7]+'_aff')
        move_mask = join(input_mask1_dir, mask3_filenames[i])
        # mov_bsp = join(out_aff, 'result.0.nrrd')
        # out_bsp = join(tmp_time2_dir, filename[:-5]+'_bsp')
        AffBsp(mov_aff, out_aff, move_mask)

        shutil.move(join(out_aff, 'result.0.nii'), join(output_time1_dir, time3_filenames[i]))

        # time3-->time1
        mov_aff = join(input_time2_dir, time3_filenames[i])
        out_aff = join(tmp_time2_dir, time3_filenames[i][:-7]+'_aff')
        move_mask = join(input_mask2_dir, mask3_filenames[i])
        # mov_bsp = join(out_aff, 'result.0.nrrd')
        # out_bsp = join(tmp_time3_dir, filename[:-5]+'_bsp')
        AffBsp(mov_aff, out_aff, move_mask)

        shutil.move(join(out_aff, 'result.0.nii'), join(output_time2_dir, time3_filenames[i]))

