import SimpleITK as sitk
import os

def pathExist(path):
    if not os.path.exists(path):
        os.makedirs(path)

    return path

path = r"F:\Zhejiang center\original image\"
out = r"F:\Zhejiang center\substraction image\"


for idx, patient in enumerate(os.listdir(path)):
    print(idx, patient)
    patient_dir = os.path.join(path, patient)
    # for idxx, img in enumerate(os.listdir(patient_dir)):
    #     print(idxx, img)
    a_file = patient_dir + '/' + 'a_img.nii.gz'
    d_file = patient_dir + '/' + 'd_img.nii.gz'
    n_file = patient_dir + '/' + 'n_img.nii.gz'
    v_file = patient_dir + '/' + 'v_img.nii.gz'

    # 读取数据
    a_img = sitk.ReadImage(a_file)
    d_img = sitk.ReadImage(d_file)
    n_img = sitk.ReadImage(n_file)
    v_img = sitk.ReadImage(v_file)

    a_data = sitk.GetArrayFromImage(a_img)
    d_data = sitk.GetArrayFromImage(d_img)
    n_data = sitk.GetArrayFromImage(n_img)
    v_data = sitk.GetArrayFromImage(v_img)

    # 图像相减
    an = a_data - n_data
    vn = v_data - n_data
    dn = d_data - n_data
    va = v_data - a_data
    da = d_data - a_data
    dv = d_data - v_data

    # 保存为减影图像
    an_img = sitk.GetImageFromArray(an)
    vn_img = sitk.GetImageFromArray(vn)
    dn_img = sitk.GetImageFromArray(dn)
    va_img = sitk.GetImageFromArray(va)
    da_img = sitk.GetImageFromArray(da)
    dv_img = sitk.GetImageFromArray(dv)

    an_img.CopyInformation(a_img)
    vn_img.CopyInformation(v_img)
    dn_img.CopyInformation(d_img)
    va_img.CopyInformation(v_img)
    da_img.CopyInformation(d_img)
    dv_img.CopyInformation(d_img)

    sitk.WriteImage(an_img, pathExist(os.path.join(out, patient)) + '/' + 'an_img.nii.gz')
    sitk.WriteImage(vn_img, pathExist(os.path.join(out, patient)) + '/' + 'vn_img.nii.gz')
    sitk.WriteImage(dn_img, pathExist(os.path.join(out, patient)) + '/' + 'dn_img.nii.gz')
    sitk.WriteImage(va_img, pathExist(os.path.join(out, patient)) + '/' + 'va_img.nii.gz')
    sitk.WriteImage(da_img, pathExist(os.path.join(out, patient)) + '/' + 'da_img.nii.gz')
    sitk.WriteImage(dv_img, pathExist(os.path.join(out, patient)) + '/' + 'dv_img.nii.gz')


# file1_path = ''
# file2_path = ''
# store_path = ''
#
# image1 = sitk.ReadImage(file1_path)
# image2 = sitk.ReadImage(file2_path)
# data1 = sitk.GetArrayFromImage(image1)
# data2 = sitk.GetArrayFromImage(image2)
# delta_data = data1 - data2
# delta_image = sitk.GetImageFromArray(delta_data)
# delta_image.CopyInformation(image1)
# sitk.WriteImage(delta_image, store_path)