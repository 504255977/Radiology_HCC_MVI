import cv2
import nibabel as nib
import os
import numpy as np
import math

def pathExist(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def dilate5():   #(0.5, 0.5, 1.0)
    image_path = r"F:\Zhejiang center\tumor mask"
    out_path   = r"F:\Zhejiang center\tumor_5mm mask"
    patient_list = os.listdir(image_path)
    for idx, patient in enumerate(patient_list):
        if idx >= 0:
            print(idx, "===", patient)
            nii_list = [i for i in os.listdir(os.path.join(image_path, patient)) if 'v_roi' in i]
            for idxxx, nii in enumerate(nii_list):
                img_nii = nib.load(os.path.join(image_path, patient, nii))
                img_data = img_nii.get_fdata()
                affine = img_nii.affine
                n = abs(round(affine[0][0], 3))
                # print(abs(round(affine[0][0], 3)))
                m = math.floor(5 * (1/n))                    # 5mm dilation of the tumor boundaries
                new_data = np.zeros(img_data.shape)
                for idx in range(min(new_data.shape)):
                    # print(idx)
                    lay = img_data[..., idx]
                    kernel = np.ones((3, 3), np.uint8)

                    dst = cv2.dilate(lay, kernel, iterations=m)  # removing some small areas

                    new_data[..., idx] = dst

                out_path1 = os.path.join(pathExist(os.path.join(out_path, patient)), nii)
                nib.save(nib.Nifti1Image(new_data, affine), out_path1)
                
def sub_mask():
    mask2 = r"F:\Zhejiang center\tumor mask"
    mask3 = r"F:\Zhejiang center\tumor_5mm mask"
    out_path = r"F:\Zhejiang center\5mm submask"
    patient_list = os.listdir(mask2)
    for idx, patient in enumerate(patient_list):
        if idx >= 0:
            print(idx, "===", patient)
            nii_list = [i for i in os.listdir(os.path.join(mask2, patient)) if 'v_roi' in i]
            for idxxx, nii in enumerate(nii_list):
                img_nii_2 = nib.load(os.path.join(mask2, patient, nii))
                img_data_2 = img_nii_2.get_fdata()
                affine = img_nii_2.affine



                img_nii_3 = nib.load(os.path.join(mask3, patient, nii))          # .replace('seg','seg_dil')))
                img_data_3 = img_nii_3.get_fdata()

                new_data = np.zeros(img_data_2.shape)

                for idx in range(min(img_data_2.shape)):
                    # print(idx)
                    dst = img_data_3[..., idx] - img_data_2[..., idx]
                    dst[dst != 1] = 0

                    new_data[..., idx] = dst
                out_path1 = os.path.join(out_path, patient)
                pathExist(out_path1)
                nib.save(nib.Nifti1Image(new_data, affine), os.path.join(out_path1, nii.split('.')[0] + '_r12mm_b6.nii.gz')) ########

if __name__ == '__main__':
    dilate5()
    #submask()
