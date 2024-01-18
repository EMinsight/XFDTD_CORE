import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import os

def graph_check(a,b, title:str = ''):
  print(title)
  print(f"Shape: {a.shape} {b.shape}")
  print(f"Max: {np.max(np.abs(a-b))}")
  max_idx = np.argmax(np.abs(a-b))
  print(f"Max a relative err: np.abs(a-b)[{max_idx}]/np.abs(a)[{max_idx}] = {np.abs(a-b)[max_idx]/np.abs(a)[max_idx]}")
  print(f"Max b relative err: np.abs(a-b)[{max_idx}]/np.abs(b)[{max_idx}] = {np.abs(a-b)[max_idx]/np.abs(b)[max_idx]}")
  print(f"Max relative : {np.max(np.abs(a-b) / np.abs(a))}")
  print(f"Min: {np.min(np.abs(a-b))}")
  print(f"Mean: {np.mean(np.abs(a-b))}")
  print(f"a unique: {np.unique(np.abs(a))}")
  print(f"b unique: {np.unique(np.abs(b))}")
  print(f"a mean: {np.mean(np.abs(a))}")
  print(f"b mean: {np.mean(np.abs(b))}")
  x = input()
  # plt.plot(a-b)
  # plt.title(title)
  # plt.show()

def check_coe(data_dir:str):
  cexe = np.load(os.path.join(data_dir,'cexe.npy'))
  cexhy = np.load(os.path.join(data_dir,'cexhy.npy'))
  cexhz = np.load(os.path.join(data_dir,'cexhz.npy'))
  ceye = np.load(os.path.join(data_dir,'ceye.npy'))
  ceyhx = np.load(os.path.join(data_dir,'ceyhx.npy'))
  ceyhz = np.load(os.path.join(data_dir,'ceyhz.npy'))
  ceze = np.load(os.path.join(data_dir,'ceze.npy'))
  cezhx = np.load(os.path.join(data_dir,'cezhx.npy'))
  cezhy = np.load(os.path.join(data_dir,'cezhy.npy'))
  chxh = np.load(os.path.join(data_dir,'chxh.npy'))
  chxey = np.load(os.path.join(data_dir,'chxey.npy'))
  chxez = np.load(os.path.join(data_dir,'chxez.npy'))
  chyh = np.load(os.path.join(data_dir,'chyh.npy'))
  chyex = np.load(os.path.join(data_dir,'chyex.npy'))
  chyez = np.load(os.path.join(data_dir,'chyez.npy'))
  chzh = np.load(os.path.join(data_dir,'chzh.npy'))
  chzex = np.load(os.path.join(data_dir,'chzex.npy'))
  chzey = np.load(os.path.join(data_dir,'chzey.npy'))

  m_cexe = sio.loadmat(os.path.join(data_dir,'m_cexe.mat'))['Cexe']
  m_cexhy = sio.loadmat(os.path.join(data_dir,'m_cexhy.mat'))['Cexhy']
  m_cexhz = sio.loadmat(os.path.join(data_dir,'m_cexhz.mat'))['Cexhz']
  m_ceye = sio.loadmat(os.path.join(data_dir,'m_ceye.mat'))['Ceye']
  m_ceyhx = sio.loadmat(os.path.join(data_dir,'m_ceyhx.mat'))['Ceyhx']
  m_ceyhz = sio.loadmat(os.path.join(data_dir,'m_ceyhz.mat'))['Ceyhz']
  m_ceze = sio.loadmat(os.path.join(data_dir,'m_ceze.mat'))['Ceze']
  m_cezhx = sio.loadmat(os.path.join(data_dir,'m_cezhx.mat'))['Cezhx']
  m_cezhy = sio.loadmat(os.path.join(data_dir,'m_cezhy.mat'))['Cezhy']
  m_chxh = sio.loadmat(os.path.join(data_dir,'m_chxh.mat'))['Chxh']
  m_chxey = sio.loadmat(os.path.join(data_dir,'m_chxey.mat'))['Chxey']
  m_chxez = sio.loadmat(os.path.join(data_dir,'m_chxez.mat'))['Chxez']
  m_chyh = sio.loadmat(os.path.join(data_dir,'m_chyh.mat'))['Chyh']
  m_chyex = sio.loadmat(os.path.join(data_dir,'m_chyex.mat'))['Chyex']
  m_chyez = sio.loadmat(os.path.join(data_dir,'m_chyez.mat'))['Chyez']
  m_chzh = sio.loadmat(os.path.join(data_dir,'m_chzh.mat'))['Chzh']
  m_chzex = sio.loadmat(os.path.join(data_dir,'m_chzex.mat'))['Chzex']
  m_chzey = sio.loadmat(os.path.join(data_dir,'m_chzey.mat'))['Chzey']

  # eps_x = np.load(os.path.join(data_dir,'eps_x.npy'))
  # m_eps_x = sio.loadmat(os.path.join(data_dir,'m_eps_x.mat'))['ans']
  # sigma_e_x = np.load(os.path.join(data_dir,'sigma_e_x.npy'))
  # m_sigma_e_x = sio.loadmat(os.path.join(data_dir,'m_sigma_e_x.mat'))['sigma_e_x']
  # print(eps_x.shape)
  # print(m_eps_x.shape)
  # print(sigma_e_x.shape)
  # print(m_sigma_e_x.shape)

  # graph_check(eps_x.flatten(),m_eps_x.flatten(), 'eps_x')
  # graph_check(sigma_e_x.flatten(),m_sigma_e_x.flatten(), 'sigma_e_x')
  # # exit(0)


  graph_check(cexe.flatten(),m_cexe.flatten(), 'Cexe')
  graph_check(cexhy.flatten(),m_cexhy.flatten(), 'Cexhy')
  graph_check(cexhz.flatten(),m_cexhz.flatten(), 'Cexhz')
  graph_check(ceye.flatten(),m_ceye.flatten(), 'Ceye')
  graph_check(ceyhx.flatten(),m_ceyhx.flatten(), 'Ceyhx')
  graph_check(ceyhz.flatten(),m_ceyhz.flatten(), 'Ceyhz')
  graph_check(ceze.flatten(),m_ceze.flatten(), 'Ceze')
  graph_check(cezhx.flatten(),m_cezhx.flatten(), 'Cezhx')
  graph_check(cezhy.flatten(),m_cezhy.flatten(), 'Cezhy')
  graph_check(chxh.flatten(),m_chxh.flatten(), 'Chxh')
  graph_check(chxey.flatten(),m_chxey.flatten(), 'Chxey')
  graph_check(chxez.flatten(),m_chxez.flatten(), 'Chxez')
  graph_check(chyh.flatten(),m_chyh.flatten(), 'Chyh')
  graph_check(chyex.flatten(),m_chyex.flatten(), 'Chyex')
  graph_check(chyez.flatten(),m_chyez.flatten(), 'Chyez')
  graph_check(chzh.flatten(),m_chzh.flatten(), 'Chzh')
  graph_check(chzex.flatten(),m_chzex.flatten(), 'Chzex')
  graph_check(chzey.flatten(),m_chzey.flatten(), 'Chzey')

def check_nffft(data_dir:str):
  jx_yp = np.load(os.path.join(data_dir,'jx_yp.npy'))
  jx_yn = np.load(os.path.join(data_dir,'jx_yn.npy'))
  jx_zp = np.load(os.path.join(data_dir,'jx_zp.npy'))
  jx_zn = np.load(os.path.join(data_dir,'jx_zn.npy'))
  jy_xp = np.load(os.path.join(data_dir,'jy_xp.npy'))
  jy_xn = np.load(os.path.join(data_dir,'jy_xn.npy'))
  jy_zp = np.load(os.path.join(data_dir,'jy_zp.npy'))
  jy_zn = np.load(os.path.join(data_dir,'jy_zn.npy'))
  jz_xp = np.load(os.path.join(data_dir,'jz_xp.npy'))
  jz_xn = np.load(os.path.join(data_dir,'jz_xn.npy'))
  jz_yp = np.load(os.path.join(data_dir,'jz_yp.npy'))
  jz_yn = np.load(os.path.join(data_dir,'jz_yn.npy'))
  mx_yp = np.load(os.path.join(data_dir,'mx_yp.npy'))
  mx_yn = np.load(os.path.join(data_dir,'mx_yn.npy'))
  mx_zp = np.load(os.path.join(data_dir,'mx_zp.npy'))
  mx_zn = np.load(os.path.join(data_dir,'mx_zn.npy'))
  my_xp = np.load(os.path.join(data_dir,'my_xp.npy'))
  my_xn = np.load(os.path.join(data_dir,'my_xn.npy'))
  my_zp = np.load(os.path.join(data_dir,'my_zp.npy'))
  my_zn = np.load(os.path.join(data_dir,'my_zn.npy'))
  mz_xp = np.load(os.path.join(data_dir,'mz_xp.npy'))
  mz_xn = np.load(os.path.join(data_dir,'mz_xn.npy'))
  mz_yp = np.load(os.path.join(data_dir,'mz_yp.npy'))
  mz_yn = np.load(os.path.join(data_dir,'mz_yn.npy'))

  m_jx_yp = sio.loadmat(os.path.join(data_dir,'m_jx_yp.mat'))['cjxyp']
  m_jx_yn = sio.loadmat(os.path.join(data_dir,'m_jx_yn.mat'))['cjxyn']
  m_jx_zp = sio.loadmat(os.path.join(data_dir,'m_jx_zp.mat'))['cjxzp']
  m_jx_zn = sio.loadmat(os.path.join(data_dir,'m_jx_zn.mat'))['cjxzn']
  m_jy_xp = sio.loadmat(os.path.join(data_dir,'m_jy_xp.mat'))['cjyxp']
  m_jy_xn = sio.loadmat(os.path.join(data_dir,'m_jy_xn.mat'))['cjyxn']
  m_jy_zp = sio.loadmat(os.path.join(data_dir,'m_jy_zp.mat'))['cjyzp']
  m_jy_zn = sio.loadmat(os.path.join(data_dir,'m_jy_zn.mat'))['cjyzn']
  m_jz_xp = sio.loadmat(os.path.join(data_dir,'m_jz_xp.mat'))['cjzxp']
  m_jz_xn = sio.loadmat(os.path.join(data_dir,'m_jz_xn.mat'))['cjzxn']
  m_jz_yp = sio.loadmat(os.path.join(data_dir,'m_jz_yp.mat'))['cjzyp']
  m_jz_yn = sio.loadmat(os.path.join(data_dir,'m_jz_yn.mat'))['cjzyn']

  m_mx_yp = sio.loadmat(os.path.join(data_dir,'m_mx_yp.mat'))['cmxyp']
  m_mx_yn = sio.loadmat(os.path.join(data_dir,'m_mx_yn.mat'))['cmxyn']
  m_mx_zp = sio.loadmat(os.path.join(data_dir,'m_mx_zp.mat'))['cmxzp']
  m_mx_zn = sio.loadmat(os.path.join(data_dir,'m_mx_zn.mat'))['cmxzn']
  m_my_xp = sio.loadmat(os.path.join(data_dir,'m_my_xp.mat'))['cmyxp']
  m_my_xn = sio.loadmat(os.path.join(data_dir,'m_my_xn.mat'))['cmyxn']
  m_my_zp = sio.loadmat(os.path.join(data_dir,'m_my_zp.mat'))['cmyzp']
  m_my_zn = sio.loadmat(os.path.join(data_dir,'m_my_zn.mat'))['cmyzn']
  m_mz_xp = sio.loadmat(os.path.join(data_dir,'m_mz_xp.mat'))['cmzxp']
  m_mz_xn = sio.loadmat(os.path.join(data_dir,'m_mz_xn.mat'))['cmzxn']
  m_mz_yp = sio.loadmat(os.path.join(data_dir,'m_mz_yp.mat'))['cmzyp']
  m_mz_yn = sio.loadmat(os.path.join(data_dir,'m_mz_yn.mat'))['cmzyn']


  graph_check(jx_yp.flatten(),m_jx_yp.flatten(), 'jx_yp')
  graph_check(jx_yn.flatten(),m_jx_yn.flatten(), 'jx_yn')
  graph_check(jx_zp.flatten(),m_jx_zp.flatten(), 'jx_zp')
  graph_check(jx_zn.flatten(),m_jx_zn.flatten(), 'jx_zn')
  graph_check(jy_xp.flatten(),m_jy_xp.flatten(), 'jy_xp')
  graph_check(jy_xn.flatten(),m_jy_xn.flatten(), 'jy_xn')
  graph_check(jy_zp.flatten(),m_jy_zp.flatten(), 'jy_zp')
  graph_check(jy_zn.flatten(),m_jy_zn.flatten(), 'jy_zn')
  graph_check(jz_xp.flatten(),m_jz_xp.flatten(), 'jz_xp')
  graph_check(jz_xn.flatten(),m_jz_xn.flatten(), 'jz_xn')
  graph_check(jz_yp.flatten(),m_jz_yp.flatten(), 'jz_yp')
  graph_check(jz_yn.flatten(),m_jz_yn.flatten(), 'jz_yn')

  graph_check(mx_yp.flatten(),m_mx_yp.flatten(), 'mx_yp')
  graph_check(mx_yn.flatten(),m_mx_yn.flatten(), 'mx_yn')
  graph_check(mx_zp.flatten(),m_mx_zp.flatten(), 'mx_zp')
  graph_check(mx_zn.flatten(),m_mx_zn.flatten(), 'mx_zn')
  graph_check(my_xp.flatten(),m_my_xp.flatten(), 'my_xp')
  graph_check(my_xn.flatten(),m_my_xn.flatten(), 'my_xn')
  graph_check(my_zp.flatten(),m_my_zp.flatten(), 'my_zp')
  graph_check(my_zn.flatten(),m_my_zn.flatten(), 'my_zn')
  graph_check(mz_xp.flatten(),m_mz_xp.flatten(), 'mz_xp')
  graph_check(mz_xn.flatten(),m_mz_xn.flatten(), 'mz_xn')
  graph_check(mz_yp.flatten(),m_mz_yp.flatten(), 'mz_yp')
  graph_check(mz_yn.flatten(),m_mz_yn.flatten(), 'mz_yn')


if __name__ == '__main__':
  # check_coe('./data/check')
  check_nffft('./data/dielectric_resonator_antenna')

