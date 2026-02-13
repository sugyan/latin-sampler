import path from 'path'
import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react-swc'

// https://vite.dev/config/
export default defineConfig({
  base: process.env.GITHUB_PAGES ? '/latin-sampler/' : '/',
  plugins: [react()],
  optimizeDeps: {
    exclude: ['latin-sampler'],
  },
  server: {
    fs: {
      allow: [path.resolve(__dirname, '..')],
    },
  },
})
