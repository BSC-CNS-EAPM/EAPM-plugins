// types.d.ts or global.d.ts
interface NodeModule {
  hot?: {
    accept(path?: string, callback?: () => void): void;
    accept(paths: string[], callback: () => void): void;
    dispose(callback: (data: unknown) => void): void;
  };
}
